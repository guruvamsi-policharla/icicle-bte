use ark_ec::short_weierstrass::{Affine as ArkAffine, Projective as ArkProjective, SWCurveConfig};
use ark_ff::{BigInteger, PrimeField};

use icicle_core::{
    curve::{Affine as IcicleAffine, Curve, Projective as IcicleProjective},
    traits::{FieldImpl, MontgomeryConvertible},
};
use icicle_runtime::{
    memory::{DeviceVec, HostSlice},
    stream::IcicleStream,
};
use rayon::prelude::*;

//============================================================================================//
//========================= Convert single field element ark<->ICICLE ========================//
//============================================================================================//
pub fn from_ark<T, I>(ark: &T) -> I
where
    T: ark_ff::Field,
    I: FieldImpl,
{
    let mut ark_bytes = vec![];
    for base_elem in ark.to_base_prime_field_elements() {
        ark_bytes.extend_from_slice(&base_elem.into_bigint().to_bytes_le());
    }
    I::from_bytes_le(&ark_bytes)
}

pub fn to_ark<T, I>(icicle: &I) -> T
where
    T: ark_ff::Field,
    I: FieldImpl,
{
    T::from_random_bytes(&icicle.to_bytes_le()).unwrap()
}

//============================================================================================//
//============================ Transmute or copy ark scalars =================================//
//============================================================================================//

// Generic function to transmute Arkworks field elements to Icicle format and return a mutable slice
pub fn transmute_ark_to_icicle_scalars<T, I>(ark_scalars: &mut [T]) -> &mut [I]
where
    T: PrimeField,
    I: FieldImpl + MontgomeryConvertible,
{
    // SAFETY: Reinterpreting Arkworks field elements as Icicle-specific scalars
    let icicle_scalars = unsafe { &mut *(ark_scalars as *mut _ as *mut [I]) };

    let icicle_host_slice = HostSlice::from_mut_slice(&mut icicle_scalars[..]);

    // Convert from Montgomery representation using the Icicle type's conversion method
    I::from_mont(icicle_host_slice, &IcicleStream::default());

    icicle_scalars
}

pub fn ark_to_icicle_scalars_async<T, I>(ark_scalars: &[T], stream: &IcicleStream) -> DeviceVec<I>
where
    T: PrimeField,
    I: FieldImpl + MontgomeryConvertible,
{
    // SAFETY: Reinterpreting Arkworks field elements as Icicle-specific scalars
    let icicle_scalars = unsafe { &*(ark_scalars as *const _ as *const [I]) };

    // Create a HostSlice from the mutable slice
    let icicle_host_slice = HostSlice::from_slice(&icicle_scalars[..]);

    let mut icicle_scalars =
        DeviceVec::<I>::device_malloc_async(ark_scalars.len(), &stream).unwrap();
    icicle_scalars
        .copy_from_host_async(&icicle_host_slice, &stream)
        .unwrap();

    // Convert from Montgomery representation using the Icicle type's conversion method
    I::from_mont(&mut icicle_scalars, &stream);
    icicle_scalars
}

pub fn ark_to_icicle_scalars<T, I>(ark_scalars: &[T]) -> DeviceVec<I>
where
    T: PrimeField,
    I: FieldImpl + MontgomeryConvertible,
{
    ark_to_icicle_scalars_async(ark_scalars, &IcicleStream::default()) // default stream is sync
}

// Note that you can also do the following but it's slower and we prefer the result in device memory
// fn ark_to_icicle_scalars<T, I>(ark_scalars: &[T]) -> Vec<I>
// where
//     T: PrimeField,
//     I: FieldImpl,
// {
//     ark_scalars
//         .par_iter()
//         .map(|ark| from_ark(ark))
//         .collect()
// }

// Note: can convert scalars back to Ark if need to by from_mont() or to_ark()

//============================================================================================//
//============================ Convert EC points ark<->ICICLE ================================//
//============================================================================================//
// Note: the generics are currently only defined for short weierstrass curves
pub fn ark_to_icicle_affine_points<AC: SWCurveConfig, IC: Curve>(
    ark_affine: &[ArkAffine<AC>],
) -> Vec<IcicleAffine<IC>>
where
    AC::BaseField: PrimeField,
{
    ark_affine
        .par_iter()
        .map(|ark| IcicleAffine {
            x: from_ark(&ark.x),
            y: from_ark(&ark.y),
        })
        .collect()
}

pub fn ark_to_icicle_projective_points<AC: SWCurveConfig, IC: Curve>(
    ark_projective: &[ArkProjective<AC>],
) -> Vec<IcicleProjective<IC>>
where
    AC::BaseField: PrimeField,
{
    ark_projective
        .par_iter()
        .map(|ark| {
            let proj_x = ark.x * ark.z;
            let proj_z = ark.z * ark.z * ark.z;
            IcicleProjective {
                x: from_ark(&proj_x),
                y: from_ark(&ark.y),
                z: from_ark(&proj_z),
            }
        })
        .collect()
}

#[allow(unused)]
pub fn icicle_to_ark_affine_points<AC: SWCurveConfig, IC: Curve>(
    icicle_affine: &[IcicleAffine<IC>],
) -> Vec<ArkAffine<AC>>
where
    AC::BaseField: ark_ff::Field,
{
    icicle_affine
        .par_iter()
        .map(|icicle| ArkAffine::new_unchecked(to_ark(&icicle.x), to_ark(&icicle.y)))
        .collect()
}

pub fn icicle_to_ark_projective_points<AC: SWCurveConfig, IC: Curve>(
    icicle_projective: &[IcicleProjective<IC>],
) -> Vec<ArkProjective<AC>>
where
    AC::BaseField: ark_ff::Field,
{
    icicle_projective
        .par_iter()
        .map(|icicle| {
            let proj_x: AC::BaseField = to_ark(&icicle.x);
            let proj_y: AC::BaseField = to_ark(&icicle.y);
            let proj_z = to_ark(&icicle.z);

            // conversion between projective used in icicle and Jacobian used in arkworks
            let proj_x = proj_x * proj_z;
            let proj_y = proj_y * proj_z * proj_z;
            ArkProjective::new_unchecked(proj_x, proj_y, proj_z)
        })
        .collect()
}
