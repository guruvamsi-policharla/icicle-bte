use crate::icicle_utils::{icicle_to_ark_projective_points, to_ark};
use crate::utils::{add_to_transcript, hash_to_bytes, xor};
use icicle_bls12_381::curve::{
    CurveCfg, G1Affine, G1Projective as G1, G2CurveCfg, G2Projective as G2, ScalarCfg, ScalarField,
};
use icicle_core::curve::Curve;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::GenerateRandom;
use merlin::Transcript;

pub struct DLogProof {
    pub c: ScalarField,       //challenge
    pub z_alpha: ScalarField, //opening for alpha
    pub z_beta: ScalarField,  //opening for beta
    pub z_s: ScalarField,     //opening for s
}

pub struct Ciphertext {
    pub ct1: [u8; 32],
    pub ct2: G2,
    pub ct3: G2,
    pub ct4: G2,
    pub gs: G1,
    pub x: ScalarField,
    pub pi: DLogProof,
}

impl Ciphertext {
    /// panicks if ciphertext does not verify
    pub fn verify(&self, htau: G2, pk: G2, g: G1, h: G2) {
        // k2.ct2^c = h^{(tau-x)*z_alpha}, k3.ct3^c = h^{z_alpha} * pk^{z_beta}, k4.ct4^c = h^{z_beta}, and k_s.gs^c = g^{z_s}
        let minus_c = ScalarField::zero() - self.pi.c;
        let recovered_k2 = (htau - (h * self.x)) * self.pi.z_alpha + (self.ct2 * minus_c);
        let recovered_k3 = h * self.pi.z_alpha + pk * self.pi.z_beta + (self.ct3 * minus_c);
        let recovered_k4 = h * self.pi.z_beta + (self.ct4 * minus_c);
        let recovered_k_s = g * self.pi.z_s + (self.gs * minus_c);

        let mut ts: Transcript = Transcript::new(&[0u8]);
        add_to_transcript(&mut ts, b"ct1", self.ct1);

        add_to_transcript(
            &mut ts,
            b"ct2",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[self.ct2])
                [0],
        );
        add_to_transcript(
            &mut ts,
            b"ct3",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[self.ct3])
                [0],
        );
        add_to_transcript(
            &mut ts,
            b"ct4",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[self.ct4])
                [0],
        );
        add_to_transcript(
            &mut ts,
            b"gs",
            icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[self.gs])[0],
        );
        add_to_transcript(
            &mut ts,
            b"x",
            to_ark::<ark_bls12_381::Fr, ScalarField>(&self.x),
        );

        add_to_transcript(
            &mut ts,
            b"k2",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[
                recovered_k2,
            ])[0],
        );
        add_to_transcript(
            &mut ts,
            b"k3",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[
                recovered_k3,
            ])[0],
        );
        add_to_transcript(
            &mut ts,
            b"k4",
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[
                recovered_k4,
            ])[0],
        );
        add_to_transcript(
            &mut ts,
            b"k_s",
            icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[
                recovered_k_s,
            ])[0],
        );

        // Fiat-Shamir to get challenge
        let mut c_bytes = [0u8; 31];
        ts.challenge_bytes(&[8u8], &mut c_bytes);
        let c = ScalarField::from_bytes_le(&c_bytes);

        // assert that the recomputed challenge matches
        assert_eq!(self.pi.c, c);
    }
}

pub fn encrypt(
    msg: [u8; 32],
    x: ScalarField,
    hid: G1,
    htau: G2,
    pk: G2,
    g: G1,
    h: G2,
) -> Ciphertext {
    let s = ScalarCfg::generate_random(1)[0];
    let gs = g * s;
    let ark_gs = icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[gs])[0];
    let hgs = hash_to_bytes(ark_gs);
    let tg = ScalarField::from_bytes_le(&hgs);

    // compute mask
    let alpha = ScalarCfg::generate_random(1)[0];
    let beta = ScalarCfg::generate_random(1)[0];

    let lhs = (hid - (g * tg)) * alpha;
    let rhs = h;

    let mask = <ark_ec::bls12::Bls12<ark_bls12_381::Config> as ark_ec::pairing::Pairing>::pairing(
        icicle_to_ark_projective_points(&[lhs])[0],
        icicle_to_ark_projective_points(&[rhs])[0],
    ); //e(H(id)/g^tg, h)^alpha

    // let mask = E::pairing(hid - (g * tg), h) * alpha; //e(H(id)/g^tg, h)^alpha
    let hmask = hash_to_bytes(mask);

    // xor msg and hmask
    let ct1: [u8; 32] = xor(&msg, &hmask).as_slice().try_into().unwrap();
    let ct2 = (htau - (h * x)) * alpha; //h^{(tau-x)*alpha}
    let ct3 = h * alpha + pk * beta; //h^alpha * pk^beta
    let ct4 = h * beta; //h^beta

    // prove knowledge of alpha, beta, and s such that ct2  = h^{(tau-x)*alpha}, ct3 = h^alpha * pk^beta, ct4 = h^beta, and gs = g^s
    // prover sends k2 = h^{(tau-x)*r_alpha}, k3 = h^{r_alpha} * pk^{r_beta}, k4 = h^{r_beta}, and k_s = g^{r_s}
    // verifier sends a random challenge c
    // prover sends z_alpha = r_alpha + c*alpha, z_beta = r_beta + c*beta, and z_s = r_s + c*s
    // verifier checks that k2.ct2^c = h^{(tau-x)*z_alpha}, k3.ct3^c = h^{z_alpha} * pk^{z_beta}, k4.ct4^c = h^{z_beta}, and k_s.gs^c = g^{z_s}

    let r_alpha = ScalarCfg::generate_random(1)[0];
    let r_beta = ScalarCfg::generate_random(1)[0];
    let r_s = ScalarCfg::generate_random(1)[0];

    let k2 = (htau - (h * x)) * r_alpha;
    let k3 = h * r_alpha + pk * r_beta;
    let k4 = h * r_beta;
    let k_s = g * r_s;

    let mut ts: Transcript = Transcript::new(&[0u8]);
    add_to_transcript(&mut ts, b"ct1", ct1);
    add_to_transcript(
        &mut ts,
        b"ct2",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[ct2])[0],
    );
    add_to_transcript(
        &mut ts,
        b"ct3",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[ct3])[0],
    );
    add_to_transcript(
        &mut ts,
        b"ct4",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[ct4])[0],
    );
    add_to_transcript(
        &mut ts,
        b"gs",
        icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[gs])[0],
    );
    add_to_transcript(&mut ts, b"x", to_ark::<ark_bls12_381::Fr, ScalarField>(&x));

    add_to_transcript(
        &mut ts,
        b"k2",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[k2])[0],
    );
    add_to_transcript(
        &mut ts,
        b"k3",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[k3])[0],
    );
    add_to_transcript(
        &mut ts,
        b"k4",
        icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[k4])[0],
    );
    add_to_transcript(
        &mut ts,
        b"k_s",
        icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[k_s])[0],
    );

    // Fiat-Shamir to get challenge
    let mut c_bytes = [0u8; 31];
    ts.challenge_bytes(&[8u8], &mut c_bytes);
    let c = ScalarField::from_bytes_le(&c_bytes);

    let z_alpha = r_alpha + c * alpha;
    let z_beta = r_beta + c * beta;
    let z_s = r_s + c * s;

    let pi = DLogProof {
        c,
        z_alpha,
        z_beta,
        z_s,
    };

    Ciphertext {
        ct1,
        ct2,
        ct3,
        ct4,
        gs,
        x,
        pi,
    }
}

#[cfg(test)]
mod tests {
    use icicle_core::ntt::{self, initialize_domain};

    use super::*;
    use crate::dealer::Dealer;

    #[test]
    fn test_encryption() {
        let batch_size = 1 << 5;
        let n = 1 << 4;

        initialize_domain(
            ntt::get_root_of_unity::<ScalarField>((2 * batch_size).try_into().unwrap()),
            &ntt::NTTInitDomainConfig::default(),
        )
        .unwrap();

        let mut dealer = Dealer::new(batch_size, n, n / 2 - 1);
        let (crs, _) = dealer.setup();
        let pk = crs.h * dealer.sk;

        let msg = [1u8; 32];
        let rou = ntt::get_root_of_unity::<ScalarField>(batch_size.try_into().unwrap());

        let hid = CurveCfg::generate_random_projective_points(1)[0];
        let ct = encrypt(msg, rou, hid, crs.htau, pk, crs.g, crs.h);

        ct.verify(crs.htau, pk, crs.g, crs.h);
    }
}
