use std::iter;
// use crate::encryption::Ciphertext;
// use crate::{icicle_utils::icicle_to_ark_affine_points, utils::lagrange_poly};
use icicle_bls12_381::curve::{
    CurveCfg, G1Affine, G1Projective as G1, G2CurveCfg, G2Projective as G2, ScalarCfg, ScalarField,
};
use icicle_core::curve::Curve;
use icicle_core::traits::FieldImpl;
use icicle_core::{ntt, traits::GenerateRandom};
use icicle_runtime::memory::{DeviceVec, HostSlice};

use crate::utils::lagrange_interp_eval_scalar;

#[derive(Clone)]
pub struct CRS {
    pub g: G1,
    pub h: G2,

    pub powers_of_g: Vec<G1Affine>,
    pub htau: G2,

    pub y: Vec<G1>, // Preprocessed Toeplitz matrix to compute opening proofs at all points
}

/// Dealer sets up the CRS and secret shares sk. Assumes the shares are over (1..n) and the secret key is stored at 0
#[derive(Clone)]
pub struct Dealer {
    pub batch_size: usize,
    pub n: usize,
    pub t: usize, // t+1 parties need to agree to decrypt
    pub sk: ScalarField,
}

impl Dealer {
    pub fn new(batch_size: usize, n: usize, t: usize) -> Self {
        let sk = ScalarCfg::generate_random(1)[0];
        Self {
            batch_size,
            n,
            t,
            sk,
        }
    }

    pub fn setup(&mut self) -> (CRS, Vec<ScalarField>) {
        // Sample tau and compute its powers ==========================================================
        let tau = ScalarCfg::generate_random(1)[0];

        let powers_of_tau: Vec<ScalarField> =
            iter::successors(Some(ScalarField::one()), |p| Some(*p * tau))
                .take(self.batch_size)
                .collect();

        // Generators
        let g = CurveCfg::generate_random_projective_points(1)[0];
        let h = G2CurveCfg::generate_random_projective_points(1)[0];

        // todo: make this a vectorized operation
        let mut powers_of_g = vec![g; self.batch_size];
        for i in 1..self.batch_size {
            powers_of_g[i] = powers_of_g[i] * powers_of_tau[i];
        }

        let htau = h * tau;

        // Compute the Toeplitz matrix preprocessing ==================================================
        let mut top_tau = powers_of_tau.clone();
        top_tau.truncate(self.batch_size);
        top_tau.reverse();
        top_tau.resize(2 * self.batch_size, ScalarField::zero());

        let mut ntt_result = DeviceVec::<ScalarField>::device_malloc(2 * self.batch_size).unwrap();
        let ntt_cfg = ntt::NTTConfig::<ScalarField>::default();
        ntt::ntt(
            HostSlice::from_slice(&top_tau),
            ntt::NTTDir::kForward,
            &ntt_cfg,
            &mut ntt_result[..],
        )
        .unwrap();

        let mut fft_top_tau = vec![ScalarField::zero(); 2 * self.batch_size];
        ntt_result
            .copy_to_host(HostSlice::from_mut_slice(fft_top_tau.as_mut_slice()))
            .unwrap();

        // todo: vectorize and do on device
        // Compute powers of top_tau
        let mut y = vec![g; 2 * self.batch_size];
        for i in 0..2 * self.batch_size {
            y[i] = y[i] * fft_top_tau[i];
        }

        let mut sk_poly = vec![ScalarField::zero(); self.t + 1];
        sk_poly[0] = self.sk;
        for i in 1..self.t {
            sk_poly[i] = ScalarCfg::generate_random(1)[0];
        }

        let share_domain = (1..=self.n)
            .map(|i| ScalarField::from_u32(i as u32))
            .collect::<Vec<_>>();

        let eval_domain = (0..=self.t)
            .map(|i| ScalarField::zero() - ScalarField::from_u32(i as u32))
            .collect::<Vec<_>>();

        let sk_shares = lagrange_interp_eval_scalar(&eval_domain, &share_domain, &sk_poly);

        // let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        // share_domain.fft_in_place(&mut sk_shares);

        let powers_of_g = powers_of_g.iter().map(|&g| g.into()).collect();

        let crs = CRS {
            g,
            h,
            powers_of_g,
            htau,
            y,
        };

        (crs, sk_shares)
    }
}

#[cfg(test)]
mod tests {
    use icicle_core::ntt::initialize_domain;

    use crate::utils::lagrange_interp_eval_g2;

    use super::*;

    #[test]
    fn test_dealer() {
        let batch_size = 1 << 5;
        let n = 1 << 4;
        let t = n / 2 - 1;

        initialize_domain(
            ntt::get_root_of_unity::<ScalarField>((2 * batch_size).try_into().unwrap()),
            &ntt::NTTInitDomainConfig::default(),
        )
        .unwrap();

        let mut dealer = Dealer::new(batch_size, n, t);
        let (crs, sk_shares) = dealer.setup();

        let share_domain = (1..=n)
            .map(|i| ScalarField::from_u32(i as u32))
            .collect::<Vec<_>>();
        let should_be_sk =
            lagrange_interp_eval_scalar(&share_domain, &vec![ScalarField::zero()], &sk_shares)[0];
        assert_eq!(dealer.sk, should_be_sk);

        let pk = crs.h * dealer.sk;
        let should_be_pk = crs.h * should_be_sk;
        assert_eq!(pk, should_be_pk);

        let g_sk_shares = sk_shares.iter().map(|&ski| crs.h * ski).collect::<Vec<_>>();

        let interp_pk =
            lagrange_interp_eval_g2(&share_domain, &vec![ScalarField::zero()], &g_sk_shares)[0];
        assert_eq!(pk, interp_pk);

        assert_eq!(crs.powers_of_g.len(), batch_size);
        assert_eq!(sk_shares.len(), n);
    }
}
