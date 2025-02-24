use crate::icicle_utils::icicle_to_ark_projective_points;
use crate::utils::commit_poly_evals;
use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, lagrange_interp_eval_g1, open_all_values, xor},
};
use ark_std::{end_timer, start_timer};
use icicle_bls12_381::curve::{
    CurveCfg, G1Projective as G1, G2CurveCfg, G2Projective as G2, ScalarField,
};
use icicle_core::ntt;
use icicle_core::traits::FieldImpl;
use icicle_runtime::memory::{DeviceVec, HostSlice};
use rayon::prelude::*;
use std::collections::BTreeMap;
pub struct SecretKey {
    sk_share: ScalarField,
}

impl SecretKey {
    pub fn new(sk_share: ScalarField) -> Self {
        SecretKey { sk_share }
    }

    /// each party in the committee computes a partial decryption
    pub fn partial_decrypt(&self, ct: &Vec<Ciphertext>, hid: G1, pk: G2, crs: &CRS) -> G1 {
        let batch_size = crs.powers_of_g.len();
        (0..batch_size).into_par_iter().for_each(|i| {
            ct[i].verify(crs.htau, pk, crs.g, crs.h);
        });

        let mut fevals = vec![ScalarField::zero(); batch_size];
        for i in 0..batch_size {
            let tg_bytes = hash_to_bytes(
                icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[ct[i].gs])
                    [0],
            );
            fevals[i] = ScalarField::from_bytes_le(&tg_bytes[0..31]);
        }

        let com = commit_poly_evals(&crs, &fevals);
        let delta = hid - com;

        let pd = delta * self.sk_share;

        pd
    }
}

/// aggregate partial decryptions into a signature on H(id)/com
pub fn aggregate_partial_decryptions(partial_decryptions: &BTreeMap<usize, G1>) -> G1 {
    // interpolate partial decryptions to recover the signature
    let mut evals = Vec::new();
    let mut eval_points = Vec::new();
    // Iterate over the map and collect keys and values
    for (&key, &value) in partial_decryptions.iter() {
        evals.push(value);
        eval_points.push(ScalarField::from_u32(key as u32));
    }

    let sigma = lagrange_interp_eval_g1(&eval_points, &vec![ScalarField::zero()], &evals)[0];

    sigma
}

/// decrypts all the ciphertexts in a batch
pub fn decrypt_all(sigma: G1, ct: &Vec<Ciphertext>, hid: G1, crs: &CRS) -> Vec<[u8; 32]> {
    let batch_size = ct.len();

    // compute fevals by hashing gs of the ciphertexts to get fevals
    let mut fevals = vec![ScalarField::zero(); batch_size];
    for i in 0..batch_size {
        let tg_bytes = hash_to_bytes(
            icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[ct[i].gs])[0],
        );
        fevals[i] = ScalarField::from_bytes_le(&tg_bytes[0..31]);
    }

    let com_timer = start_timer!(|| "Commit fevals");
    let com = commit_poly_evals(&crs, &fevals);
    end_timer!(com_timer);

    let delta = hid - com;

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    let fk22_timer = start_timer!(|| "FK22");
    let mut ntt_result = DeviceVec::<ScalarField>::device_malloc(batch_size).unwrap();
    let ntt_cfg = ntt::NTTConfig::<ScalarField>::default();
    ntt::ntt(
        HostSlice::from_slice(&fevals),
        ntt::NTTDir::kInverse,
        &ntt_cfg,
        &mut ntt_result[..],
    )
    .unwrap();

    let mut fcoeffs = vec![ScalarField::zero(); batch_size];
    ntt_result
        .copy_to_host(HostSlice::from_mut_slice(fcoeffs.as_mut_slice()))
        .unwrap();

    let pi = open_all_values(&crs.y, &fcoeffs, batch_size);
    end_timer!(fk22_timer);

    // now decrypt each of the ciphertexts as m = ct(1) xor H(e(delta,ct2).e(pi,ct3)e(-sigma,ct4))
    let ppe_timer = start_timer!(|| "PPEs");
    let mut m = vec![[0u8; 32]; batch_size];

    let ark_pi = icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&pi);
    let ark_delta =
        icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[delta])[0];
    let ark_sigma =
        icicle_to_ark_projective_points::<ark_bls12_381::g1::Config, CurveCfg>(&[sigma])[0];

    let ark_ct = ct
        .par_iter()
        .map(|c| {
            icicle_to_ark_projective_points::<ark_bls12_381::g2::Config, G2CurveCfg>(&[
                c.ct2, c.ct3, c.ct4,
            ])
        })
        .collect::<Vec<_>>();

    m.par_iter_mut().enumerate().for_each(|(i, m_i)| {
        let ark_g1_terms = vec![ark_pi[i], ark_delta, - ark_sigma];
        let ark_g2_terms = vec![ark_ct[i][0], ark_ct[i][1], ark_ct[i][2]];
            
        let mask = <ark_ec::bls12::Bls12<ark_bls12_381::Config> as ark_ec::pairing::Pairing>::multi_miller_loop(
            &ark_g1_terms,
            &ark_g2_terms,
        );
        let mask = <ark_ec::bls12::Bls12<ark_bls12_381::Config> as ark_ec::pairing::Pairing>::final_exponentiation(mask).unwrap();

        let hmask = hash_to_bytes(mask);
        *m_i = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    });
    end_timer!(ppe_timer);

    m
}
