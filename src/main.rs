use clap::Parser;
use icicle_bls12_381::curve::{ScalarCfg, ScalarField};
use icicle_bte::{dealer::Dealer, utils::open_all_values};
use icicle_core::{
    ntt::{self, initialize_domain},
    traits::GenerateRandom,
};
use std::convert::TryInto;
use std::time::Instant;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long, default_value_t = 14)]
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

    let args = Args::parse();
    println!("{:?}", args);

    try_load_and_set_backend_device(&args);

    let max_log_size = args.size;
    let domain_size = 1 << max_log_size;
    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * domain_size).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    /////////////////////////////////

    let n = 1 << 4;
    let sample_size = 1;

    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * (1 << max_log_size)).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    for size in max_log_size..=max_log_size {
        let batch_size = 1 << size;

        let mut dealer = Dealer::new(batch_size, n, n / 2 - 1);
        let (crs, _) = dealer.setup();

        let f_coeffs = ScalarCfg::generate_random(batch_size);

        // bench full decryption
        let start = Instant::now();
        for _ in 0..sample_size {
            open_all_values(&crs.y, &f_coeffs, batch_size);
        }
        let time_taken = start.elapsed().as_millis() as f64;
        println!(
            "open_all_values: batch_size = {}, time = {}ms",
            batch_size,
            time_taken / (sample_size as f64)
        );
    }
}
