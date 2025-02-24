use ark_std::{end_timer, start_timer};
use clap::Parser;
use icicle_bls12_381::curve::{CurveCfg, G1Projective as G1, ScalarField};
use icicle_bte::{
    dealer::Dealer,
    decryption::{aggregate_partial_decryptions, decrypt_all, SecretKey},
    encryption::{encrypt, Ciphertext},
};
use icicle_core::ntt::{self, initialize_domain};
use icicle_core::{curve::Curve, traits::FieldImpl};
use rayon::prelude::*;
use std::collections::BTreeMap;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long, default_value_t = 10)]
    size: u8,

    /// Device type (e.g., "CPU", "CUDA")
    #[arg(short, long, default_value = "CPU")]
    device_type: String,
}

// Load backend and set device
fn try_load_and_set_backend_device(args: &Args) {
    if args.device_type != "CPU" {
        icicle_runtime::runtime::load_backend_from_env_or_default().unwrap();
    }
    println!("Setting device {}", args.device_type);
    let device = icicle_runtime::Device::new(&args.device_type, 0 /* =device_id*/);
    icicle_runtime::set_device(&device).unwrap();
}

fn main() {
    let args = Args::parse();
    println!("{:?}", args);

    try_load_and_set_backend_device(&args);

    let batch_size = 1 << args.size;
    let n = 1 << 2;

    println!("n: {}, batch_size: {}", n, batch_size);

    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * batch_size).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    println!("Generating CRS and secret keys");
    let now = std::time::Instant::now();
    let mut dealer = Dealer::new(batch_size, n, n / 2 - 1);
    let (crs, sk_shares) = dealer.setup();
    let pk = crs.h * dealer.sk;

    let mut secret_key: Vec<SecretKey> = Vec::new();
    for i in 0..n {
        secret_key.push(SecretKey::new(sk_shares[i]));
    }
    println!(
        "Total Time taken for CRS and secret key generation: {}ms",
        now.elapsed().as_millis()
    );

    let msg = [1u8; 32];
    let rou = ntt::get_root_of_unity::<ScalarField>(batch_size.try_into().unwrap());
    let mut domain_elements: Vec<ScalarField> = vec![ScalarField::one()];
    for i in 1..batch_size {
        domain_elements.push(domain_elements[i - 1] * rou);
    }

    let hid = CurveCfg::generate_random_projective_points(1)[0];

    // generate ciphertexts for all points in tx_domain
    println!("Encrypting messages");
    let ct: Vec<Ciphertext> = domain_elements
        .par_iter()
        .map(|&x| encrypt(msg, x, hid, crs.htau, pk, crs.g, crs.h))
        .collect();

    // generate partial decryptions
    println!("Generating partial decryptions");
    let mut partial_decryptions: BTreeMap<usize, G1> = BTreeMap::new();
    for i in 0..n / 2 {
        let partial_decryption = secret_key[i].partial_decrypt(&ct, hid, pk, &crs);
        partial_decryptions.insert(i + 1, partial_decryption);
    }

    // aggregate partial decryptions
    println!("Aggregating partial decryptions");
    let sigma = aggregate_partial_decryptions(&partial_decryptions);

    // verify that ciphertexts are well formed
    println!("Verifying ciphertexts");
    ct.par_iter().for_each(|ciphertext| {
        ciphertext.verify(crs.htau, pk, crs.g, crs.h);
    });

    let dec_timer = start_timer!(|| "Decryption");
    let messages = decrypt_all(sigma, &ct, hid, &crs);
    end_timer!(dec_timer);

    for i in 0..batch_size {
        assert_eq!(msg, messages[i]);
    }
}
