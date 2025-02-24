use crate::dealer::CRS;
use ark_serialize::CanonicalSerialize;
use icicle_bls12_381::curve::{G1Projective as G1, ScalarField};
use icicle_bls12_381::polynomials::DensePolynomial;
use icicle_core::ecntt::ecntt_inplace;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::traits::{Arithmetic, FieldImpl};
use icicle_core::{msm, ntt};
use icicle_runtime::memory::{DeviceVec, HostSlice};
use icicle_runtime::stream::IcicleStream;
use merlin::Transcript;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::ops::Div;

pub fn hash_to_bytes<T: CanonicalSerialize>(inp: T) -> [u8; 32] {
    let mut bytes = Vec::new();
    inp.serialize_compressed(&mut bytes).unwrap();
    let hash = blake3::hash(bytes.as_slice());
    let hash_bytes = hash.as_bytes();
    *hash_bytes
}

pub fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x ^ y).collect()
}

pub fn add_to_transcript<T: CanonicalSerialize>(
    ts: &mut Transcript,
    label: &'static [u8],
    data: T,
) {
    let mut data_bytes = Vec::new();
    data.serialize_compressed(&mut data_bytes).unwrap();
    ts.append_message(label, &data_bytes);
}

/// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
pub fn lagrange_poly(n: usize, i: usize) -> DensePolynomial {
    debug_assert!(i < n);
    //todo: check n is a power of 2
    let mut evals = vec![ScalarField::zero(); n];
    evals[i] = ScalarField::one();

    let poly = DensePolynomial::from_rou_evals(HostSlice::from_slice(&evals), n);

    poly
}

/// given evaluations of a polynomial over given_domain
/// interpolates the polynomial and evaluates it on target_domain
pub fn lagrange_interp_eval_g1(
    given_domain: &Vec<ScalarField>,
    target_domain: &Vec<ScalarField>,
    evals: &Vec<G1>,
) -> Vec<G1> {
    debug_assert_eq!(
        given_domain.len(),
        evals.len(),
        "Evals length does not match given_domain length"
    );

    let mut result = Vec::new();
    for &point in target_domain.iter() {
        let mut lagrange_coeffs = vec![ScalarField::one(); given_domain.len()];

        for i in 0..given_domain.len() {
            let mut num = ScalarField::one();
            let mut denom = ScalarField::one();
            for j in 0..given_domain.len() {
                if given_domain[i] != given_domain[j] {
                    num = num * (point - given_domain[j]);
                    denom = denom * (given_domain[i] - given_domain[j]);
                }
            }
            lagrange_coeffs[i] = num * denom.inv();
        }

        let mut point_eval = G1::zero();
        for i in 0..given_domain.len() {
            point_eval = point_eval + evals[i] * lagrange_coeffs[i];
        }

        result.push(point_eval);
    }

    result
}

/// given evaluations of a polynomial over given_domain
/// interpolates the polynomial and evaluates it on target_domain
pub fn lagrange_interp_eval_scalar(
    given_domain: &Vec<ScalarField>,
    target_domain: &Vec<ScalarField>,
    evals: &Vec<ScalarField>,
) -> Vec<ScalarField> {
    debug_assert_eq!(
        given_domain.len(),
        evals.len(),
        "Evals length does not match given_domain length"
    );

    let mut result = Vec::new();
    for &point in target_domain.iter() {
        let mut lagrange_coeffs = vec![ScalarField::one(); given_domain.len()];

        for i in 0..given_domain.len() {
            let mut num = ScalarField::one();
            let mut denom = ScalarField::one();
            for j in 0..given_domain.len() {
                if given_domain[i] != given_domain[j] {
                    num = num * (point - given_domain[j]);
                    denom = denom * (given_domain[i] - given_domain[j]);
                }
            }
            lagrange_coeffs[i] = num * denom.inv();
        }

        let mut point_eval = ScalarField::zero();
        for i in 0..given_domain.len() {
            point_eval = point_eval + evals[i] * lagrange_coeffs[i];
        }

        result.push(point_eval);
    }

    result
}

pub fn commit_poly(crs: &CRS, polynomial: &mut DensePolynomial) -> G1 {
    let mut com_on_device = DeviceVec::<G1>::device_malloc(1).unwrap();
    let deg = polynomial.degree() as usize;

    let mut g1_stream = IcicleStream::create().unwrap();
    let mut g1_cfg = msm::MSMConfig::default();
    g1_cfg.stream_handle = *g1_stream;
    g1_cfg.is_async = true;

    msm::msm(
        &polynomial.coeffs_mut_slice()[0..=deg as usize],
        HostSlice::from_slice(&crs.powers_of_g[0..=deg as usize]),
        &g1_cfg,
        &mut com_on_device[..],
    )
    .unwrap();

    let mut com = vec![G1::zero()];
    com_on_device
        .copy_to_host(HostSlice::from_mut_slice(com.as_mut_slice()))
        .unwrap();

    g1_stream.destroy().unwrap();

    com[0]
}

pub fn commit_poly_evals(crs: &CRS, evals: &Vec<ScalarField>) -> G1 {
    let mut com_on_device = DeviceVec::<G1>::device_malloc(1).unwrap();

    let mut g1_stream = IcicleStream::create().unwrap();
    let mut g1_cfg = msm::MSMConfig::default();
    g1_cfg.stream_handle = *g1_stream;
    g1_cfg.is_async = true;

    msm::msm(
        HostSlice::from_slice(&evals),
        HostSlice::from_slice(&crs.lagrange_powers_of_g[0..evals.len()]),
        &g1_cfg,
        &mut com_on_device[..],
    )
    .unwrap();

    let mut com = vec![G1::zero()];
    com_on_device
        .copy_to_host(HostSlice::from_mut_slice(com.as_mut_slice()))
        .unwrap();

    g1_stream.destroy().unwrap();

    com[0]
}

/// compute KZG opening proof
pub fn compute_opening_proof(crs: &CRS, polynomial: &DensePolynomial, point: ScalarField) -> G1 {
    let eval = polynomial.eval(&point);

    let mut numerator = polynomial.clone();
    numerator.sub_monomial_inplace(&eval, 0);
    let divisor = DensePolynomial::from_coeffs(
        HostSlice::from_slice(&vec![ScalarField::zero() - point, ScalarField::one()]),
        2,
    );
    let mut witness_polynomial = numerator.div(&divisor);

    // commit to b in g2
    let mut pi_on_device = DeviceVec::<G1>::device_malloc(1).unwrap();
    let deg = witness_polynomial.degree() as usize;

    let mut g1_stream = IcicleStream::create().unwrap();
    let mut g1_cfg = msm::MSMConfig::default();
    g1_cfg.stream_handle = *g1_stream;
    g1_cfg.is_async = true;

    msm::msm(
        &witness_polynomial.coeffs_mut_slice()[0..=deg as usize],
        HostSlice::from_slice(&crs.powers_of_g[0..=deg as usize]),
        &g1_cfg,
        &mut pi_on_device[..],
    )
    .unwrap();

    let mut pi = vec![G1::zero()];
    pi_on_device
        .copy_to_host(HostSlice::from_mut_slice(pi.as_mut_slice()))
        .unwrap();
    let pi = pi[0];

    g1_stream.destroy().unwrap();

    pi
}

/// Computes all the openings of a KZG commitment in O(n log n) time
/// See https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf
/// eprint version has a bug and hasn't been updated
pub fn open_all_values(y: &Vec<G1>, f: &Vec<ScalarField>, batch_size: usize) -> Vec<G1> {
    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    // f = {f0 ,f1, ..., fd}
    // v = {(d 0s), f1, ..., fd}
    let mut v = vec![ScalarField::zero(); batch_size + 1];
    v.append(&mut f[1..].to_vec());

    debug_assert_eq!(v.len(), 2 * batch_size);

    let mut v_dev = DeviceVec::<ScalarField>::device_malloc(2 * batch_size).unwrap();
    v_dev.copy_from_host(HostSlice::from_slice(&v)).unwrap();

    let ntt_cfg = ntt::NTTConfig::<ScalarField>::default();
    ntt::ntt_inplace(&mut v_dev, ntt::NTTDir::kForward, &ntt_cfg).unwrap();

    let mut fft_v = vec![ScalarField::zero(); 2 * batch_size];
    v_dev
        .copy_to_host(HostSlice::from_mut_slice(fft_v.as_mut_slice()))
        .unwrap();

    // h = y \odot v
    let mut h = vec![G1::zero(); 2 * batch_size];
    h.par_iter_mut()
        .zip(fft_v.par_iter())
        .zip(y.par_iter())
        .for_each(|((h_i, &fft_v_i), &y_i)| {
            *h_i = y_i * fft_v_i;
        });

    // inverse fft on h
    ecntt_inplace(
        HostSlice::from_mut_slice(&mut h),
        ntt::NTTDir::kInverse,
        &ntt_cfg,
    )
    .unwrap();

    h.truncate(batch_size);

    // fft on h to get KZG proofs
    ecntt_inplace(
        HostSlice::from_mut_slice(&mut h),
        ntt::NTTDir::kForward,
        &ntt_cfg,
    )
    .unwrap();

    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{dealer::Dealer, icicle_utils::icicle_to_ark_projective_points};
    use ark_ff::Zero;
    use icicle_bls12_381::curve::{G2Projective as G2, ScalarCfg};
    use icicle_core::ntt::{self, initialize_domain};
    use icicle_core::traits::GenerateRandom;

    #[test]
    fn kzg_open_test() {
        let domain_size = 1 << 5;
        initialize_domain(
            ntt::get_root_of_unity::<ScalarField>((2 * domain_size).try_into().unwrap()),
            &ntt::NTTInitDomainConfig::default(),
        )
        .unwrap();

        let mut dealer = Dealer::new(domain_size, 1 << 5, domain_size / 2 - 1);
        let (crs, _) = dealer.setup();

        let fcoeffs = ScalarCfg::generate_random(domain_size);
        let mut f = DensePolynomial::from_coeffs(HostSlice::from_slice(&fcoeffs), domain_size);

        let mut ntt_result = DeviceVec::<ScalarField>::device_malloc(domain_size).unwrap();
        let ntt_cfg = ntt::NTTConfig::<ScalarField>::default();
        ntt::ntt(
            HostSlice::from_slice(&fcoeffs),
            ntt::NTTDir::kForward,
            &ntt_cfg,
            &mut ntt_result[..],
        )
        .unwrap();

        let mut fevals = vec![ScalarField::zero(); domain_size];
        ntt_result
            .copy_to_host(HostSlice::from_mut_slice(fevals.as_mut_slice()))
            .unwrap();

        let com = commit_poly(&crs, &mut f);
        let evals_com = commit_poly_evals(&crs, &fevals);

        assert_eq!(com, evals_com);

        let point = ScalarCfg::generate_random(1)[0];

        let pi = compute_opening_proof(&crs, &f, point);

        let g1_terms = [com - (crs.g * f.eval(&point)), pi].to_vec();
        let g2_terms = [G2::zero() - crs.h, crs.htau - (crs.h * point)].to_vec();

        let ark_g1_terms: Vec<ark_bls12_381::G1Projective> =
            icicle_to_ark_projective_points(&g1_terms);

        let ark_g2_terms: Vec<ark_bls12_381::G2Projective> =
            icicle_to_ark_projective_points(&g2_terms);

        let should_be_zero = <ark_ec::bls12::Bls12::<ark_bls12_381::Config> as ark_ec::pairing::Pairing>::multi_pairing(
            ark_g1_terms,
            ark_g2_terms,
        );

        assert!(should_be_zero.is_zero());
    }

    #[test]
    fn open_all_test() {
        let domain_size = 1 << 10;
        initialize_domain(
            ntt::get_root_of_unity::<ScalarField>((2 * domain_size).try_into().unwrap()),
            &ntt::NTTInitDomainConfig::default(),
        )
        .unwrap();

        let mut dealer = Dealer::new(domain_size, 1 << 5, domain_size / 2 - 1);
        let (crs, _) = dealer.setup();

        let f_coeffs = ScalarCfg::generate_random(domain_size);
        let mut f = DensePolynomial::from_coeffs(HostSlice::from_slice(&f_coeffs), domain_size);

        let com = commit_poly(&crs, &mut f);

        let pi = open_all_values(&crs.y, &f_coeffs, domain_size);

        // compute domain elements = [1, rou, rou^2, ...] as powers of rou using dynamic programming
        let rou = ntt::get_root_of_unity::<ScalarField>(domain_size.try_into().unwrap());
        let mut domain_elements: Vec<ScalarField> = vec![ScalarField::one()];
        for i in 1..domain_size {
            domain_elements.push(domain_elements[i - 1] * rou);
        }

        // verify the kzg proof
        for i in 0..domain_size {
            let g1_terms = [com - (crs.g * f.eval(&domain_elements[i])), pi[i]].to_vec();
            let g2_terms = [G2::zero() - crs.h, crs.htau - (crs.h * domain_elements[i])].to_vec();

            let ark_g1_terms: Vec<ark_bls12_381::G1Projective> =
                icicle_to_ark_projective_points(&g1_terms);

            let ark_g2_terms: Vec<ark_bls12_381::G2Projective> =
                icicle_to_ark_projective_points(&g2_terms);

            let should_be_zero = <ark_ec::bls12::Bls12::<ark_bls12_381::Config> as ark_ec::pairing::Pairing>::multi_pairing(
            ark_g1_terms,
            ark_g2_terms,
        );

            assert!(should_be_zero.is_zero());
        }
    }

    #[test]
    fn lagrange_interp_eval_test() {
        let domain_size = 1 << 4;
        let domain = (0..domain_size)
            .map(|i| ScalarField::from_u32(i as u32))
            .collect::<Vec<_>>();

        let points = (0..domain_size / 2)
            .map(|i| ScalarField::from_u32((domain_size + i) as u32))
            .collect::<Vec<_>>();

        let f = ScalarCfg::generate_random(domain_size);

        let f = DensePolynomial::from_coeffs(HostSlice::from_slice(&f), domain_size);

        let evals = domain.iter().map(|&e| f.eval(&e)).collect::<Vec<_>>();

        let computed_evals = lagrange_interp_eval_scalar(&domain, &points, &evals);
        let should_be_evals = points.iter().map(|p| f.eval(p)).collect::<Vec<_>>();

        assert_eq!(computed_evals, should_be_evals);
    }
}
