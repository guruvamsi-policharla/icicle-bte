use icicle_bls12_381::curve::{CurveCfg, G1Projective as G1, ScalarField};
use icicle_bte::{
    dealer::Dealer,
    decryption::{aggregate_partial_decryptions, decrypt_all, SecretKey},
    encryption::{encrypt, Ciphertext},
    icicle_utils::icicle_to_ark_projective_points,
};
use icicle_core::ntt::{self, initialize_domain};
use icicle_core::{curve::Curve, traits::FieldImpl};
use std::collections::BTreeMap;

fn main() {
    let batch_size = 1 << 5;
    let n = 1 << 3;

    println!("n: {}, batch_size: {}", n, batch_size);

    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * batch_size).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    let mut dealer = Dealer::new(batch_size, n, n / 2 - 1);
    let (crs, sk_shares) = dealer.setup();
    let pk = crs.h * dealer.sk;

    let mut secret_key: Vec<SecretKey> = Vec::new();
    for i in 0..n {
        secret_key.push(SecretKey::new(sk_shares[i]));
    }

    let msg = [1u8; 32];
    let rou = ntt::get_root_of_unity::<ScalarField>(batch_size.try_into().unwrap());
    let mut domain_elements: Vec<ScalarField> = vec![ScalarField::one()];
    for i in 1..batch_size {
        domain_elements.push(domain_elements[i - 1] * rou);
    }

    let hid = CurveCfg::generate_random_projective_points(1)[0];

    // generate ciphertexts for all points in tx_domain
    let mut ct: Vec<Ciphertext> = Vec::new();
    for x in domain_elements {
        ct.push(encrypt(msg, x, hid, crs.htau, pk, crs.g, crs.h));
    }

    // generate partial decryptions
    let mut partial_decryptions: BTreeMap<usize, G1> = BTreeMap::new();
    for i in 0..n / 2 {
        let partial_decryption = secret_key[i].partial_decrypt(&ct, hid, pk, &crs);
        partial_decryptions.insert(i + 1, partial_decryption);
    }

    let sigma = aggregate_partial_decryptions(&partial_decryptions);

    // verify that ciphertexts are well formed
    for i in 0..batch_size {
        ct[i].verify(crs.htau, pk, crs.g, crs.h);
    }

    let messages = decrypt_all(sigma, &ct, hid, &crs);
    for i in 0..batch_size {
        assert_eq!(msg, messages[i]);
    }
}
