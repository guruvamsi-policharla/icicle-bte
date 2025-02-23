use icicle_bls12_381::curve::{ScalarCfg, ScalarField};
use icicle_bte::{dealer::Dealer, utils::open_all_values};
use icicle_core::ntt::{self, initialize_domain};
use icicle_core::traits::GenerateRandom;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn bench_open_all(c: &mut Criterion) {
    let max_log_size = 10;
    let domain_size = 1 << max_log_size;
    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * domain_size).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    /////////////////////////////////

    let n = 1 << 4;
    let mut group = c.benchmark_group("open_all");
    group.sample_size(20);

    let max_log_size = 10;
    initialize_domain(
        ntt::get_root_of_unity::<ScalarField>((2 * (1 << max_log_size)).try_into().unwrap()),
        &ntt::NTTInitDomainConfig::default(),
    )
    .unwrap();

    for size in 5..=max_log_size {
        let batch_size = 1 << size;

        let mut dealer = Dealer::new(batch_size, n, n / 2 - 1);
        let (crs, _) = dealer.setup();

        let f_coeffs = ScalarCfg::generate_random(batch_size);

        // bench full decryption
        group.bench_with_input(
            BenchmarkId::from_parameter(batch_size),
            &(crs.y, f_coeffs, batch_size),
            |b, inp| {
                b.iter(|| {
                    open_all_values(&inp.0, &inp.1, inp.2);
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_open_all);
criterion_main!(benches);
