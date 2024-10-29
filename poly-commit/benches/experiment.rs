use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
};
use ark_ec::{bls12::Bls12, pairing::Pairing};
use ark_pcs_bench_templates::*;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial as DenseUnivariatePoly;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

use ark_poly_commit::ipa_pc::InnerProductArgPC;
use ark_poly_commit::linear_codes::{LinearCodePCS, MultilinearBrakedown,MultilinearLigero};
use ark_poly_commit::hyrax::HyraxPC;
use ark_poly_commit::kzg10::KZG10;
use blake2::Blake2s256;

use ark_ff::PrimeField;

use rand_chacha::ChaCha20Rng;

const MIN_NUM_VARS: usize = 28;
const MAX_NUM_VARS: usize = 29;

fn main() {
    // experiment_brakedown();
    // experiment_hyrax();
    // experiment_lingero();
    experiment_ipa();
}

fn experiment_brakedown() {
    struct MerkleTreeParams;
    type LeafH = LeafIdentityHasher;
    type CompressH = Sha256;
    impl Config for MerkleTreeParams {
        type Leaf = Vec<u8>;

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    pub type MLE<F> = DenseMultilinearExtension<F>;
    type MTConfig = MerkleTreeParams;
    type ColHasher<F> = FieldToBytesColHasher<F, Blake2s256>;
    type Brakedown<F> = LinearCodePCS<
        MultilinearBrakedown<F, MTConfig, MLE<F>, ColHasher<F>>,
        F,
        MLE<F>,
        MTConfig,
        ColHasher<F>,
    >;

    fn rand_poly_brakedown_ml<F: PrimeField>(
        num_vars: usize,
        rng: &mut ChaCha20Rng,
    ) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::rand(num_vars, rng)
    }
    
    fn rand_point_brakedown_ml<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
        (0..num_vars).map(|_| F::rand(rng)).collect()
    }

    println!("\nBrakedown on BN254:");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        experiment::<ark_bn254::Fr, MLE<ark_bn254::Fr>, Brakedown<ark_bn254::Fr>>(
            num_vars,
            rand_poly_brakedown_ml,
            rand_point_brakedown_ml
        );
    }
}

fn experiment_hyrax() {
    type Hyrax254 = HyraxPC<ark_bn254::G1Affine, DenseMultilinearExtension<ark_bn254::Fr>>;

    fn rand_poly_hyrax<F: PrimeField>(
        num_vars: usize,
        rng: &mut ChaCha20Rng,
    ) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::rand(num_vars, rng)
    }

    fn rand_point_hyrax<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
        (0..num_vars).map(|_| F::rand(rng)).collect()
    }

    println!("\nHyrax on BN254:");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        experiment::<ark_bn254::Fr, DenseMultilinearExtension<ark_bn254::Fr>, Hyrax254>(
            num_vars,
            rand_poly_hyrax,
            rand_point_hyrax
        );
    }
}

fn experiment_lingero() {
    struct MerkleTreeParams;
    type LeafH = LeafIdentityHasher;
    type CompressH = Sha256;
    impl Config for MerkleTreeParams {
        type Leaf = Vec<u8>;

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    pub type MLE<F> = DenseMultilinearExtension<F>;
    type MTConfig = MerkleTreeParams;
    type ColHasher<F> = FieldToBytesColHasher<F, Blake2s256>;
    type Ligero<F> = LinearCodePCS<
        MultilinearLigero<F, MTConfig, MLE<F>, ColHasher<F>>,
        F,
        MLE<F>,
        MTConfig,
        ColHasher<F>,
    >;

    fn rand_poly_ligero_ml<F: PrimeField>(
        num_vars: usize,
        rng: &mut ChaCha20Rng,
    ) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::rand(num_vars, rng)
    }

    fn rand_point_ligero_ml<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
        (0..num_vars).map(|_| F::rand(rng)).collect()
    }

    println!("\nHyrax on BN254:");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        experiment::<ark_bn254::Fr, MLE<ark_bn254::Fr>, Ligero<ark_bn254::Fr>>(
            num_vars,
            rand_poly_ligero_ml,
            rand_point_ligero_ml
        );
    }
}

fn experiment_ipa() {

    type UniPoly = DenseUnivariatePoly<ark_ed_on_bls12_381::Fr>;
    // IPA_PC over the JubJub curve with Blake2s as the hash function
    #[allow(non_camel_case_types)]
    type IPA_JubJub = InnerProductArgPC<ark_ed_on_bls12_381::EdwardsAffine, Blake2s256, UniPoly>;

    fn rand_poly_ipa_pc<F: PrimeField>(degree: usize, rng: &mut ChaCha20Rng) -> DenseUnivariatePoly<F> {
        DenseUnivariatePoly::rand(degree, rng)
    }

    fn rand_point_ipa_pc<F: PrimeField>(_: usize, rng: &mut ChaCha20Rng) -> F {
        F::rand(rng)
    }

    println!("\nIPA on BLS12-381:");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        experiment::<ark_ed_on_bls12_381::Fr, UniPoly, IPA_JubJub>(
            2_usize.pow(21 as u32),
            rand_poly_ipa_pc,
            rand_point_ipa_pc
        );
    }
}

// fn experiment_kzg() {
//     #![allow(non_camel_case_types)]
//     use crate::kzg10::*;
//     use crate::*;

//     use ark_bls12_377::Bls12_377;
//     use ark_bls12_381::Bls12_381;
//     use ark_bls12_381::Fr;
//     use ark_poly::univariate::DensePolynomial as DensePoly;
//     use ark_std::test_rng;

//     type UniPoly_381 = DensePoly<<Bls12_381 as Pairing>::ScalarField>;
//     type UniPoly_377 = DensePoly<<Bls12_377 as Pairing>::ScalarField>;
//     type KZG_Bls12_381 = KZG10<Bls12_381, UniPoly_381>;
// }

// fn end_to_end_test_template<E, P>()
// where
//     E: Pairing,
//     P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
//     for<'a, 'b> &'a P: Div<&'b P, Output = P>,
// {
//     #![allow(non_camel_case_types)]
//     use ark_poly_commit::kzg10::*;
//     use ark_poly_commit::*;

//     use ark_bls12_377::Bls12_377;
//     use ark_bls12_381::Bls12_381;
//     use ark_bls12_381::Fr;
//     use ark_poly::univariate::DensePolynomial as DensePoly;
//     use ark_std::test_rng;

//     type UniPoly_381 = DensePoly<<Bls12_381 as Pairing>::ScalarField>;
//     type UniPoly_377 = DensePoly<<Bls12_377 as Pairing>::ScalarField>;
//     type KZG_Bls12_381 = KZG10<Bls12_381, UniPoly_381>;

//     let rng = &mut test_rng();

//     let mut degree = 0;
//     while degree <= 1 {
//         degree = usize::rand(rng) % 20;
//     }
//     let pp = KZG10::<E, P>::setup(degree, false, rng)?;
//     let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
//     let p = P::rand(degree, rng);
//     let hiding_bound = Some(1);
//     let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
//     let point = E::ScalarField::rand(rng);
//     let value = p.evaluate(&point);
//     let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
//     assert!(
//         KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
//         "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
//         degree,
//         p.degree(),
//         hiding_bound,
//     );
// }