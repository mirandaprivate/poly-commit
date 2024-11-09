use std::time::Instant;

// use core::num;

use crate::Error;
use ark_ec::pairing::{Pairing, PairingOutput};
    // scalar_mul::ScalarMul, AffineRepr, CurveGroup, VariableBaseMSM
// };
use ark_ff::{
    One,
    // PrimeField,
    Field,
    UniformRand,
    // Zero
};
// use ark_poly::DenseUVPolynomial;
use ark_std::{
    // format,
    marker::PhantomData, ops::Mul, rand::RngCore
};

use ark_serialize::CanonicalSerialize;

use rand::thread_rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;


mod data_structures;
pub use data_structures::*;

mod utils;
pub use utils::*;

/// SMART_PC is an implementation of modified Dory
pub struct SmartPC<E: Pairing> {
    _engine: PhantomData<E>,
}

impl<E> SmartPC<E>
where
    E: Pairing,
{
    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    ///
    pub fn setup<R: RngCore>(
        num_vars: usize,
        rng: &mut R,
    ) -> Result<UniversalParams<E>, Error> {
        if num_vars < 1 {
            return Err(Error::DegreeIsZero);
        }

        let q = 2_usize.pow((num_vars/2) as u32);
        let setup_time =
        start_timer!(||
            format!("Smart-PC::Setup with dimension 2^{}, 2^{} coefficients",
            num_vars,
            2 * num_vars)
        );
        
        let nu = E::ScalarField::rand(rng);
        let hat_s = E::ScalarField::rand(rng);
        let g_hat = E::G1::rand(rng);
        let h_hat = E::G2::rand(rng);


        // powers_of_s = [1, s, ..., s^(q-1)], len = q
        let mut powers_of_s = vec![E::ScalarField::one()];
        let mut cur = hat_s;
        for _ in 0..(q-1) {
            powers_of_s.push(cur);
            cur *= &hat_s;
        }

        let hat_sq = cur;

        let mut powers_of_sq = vec![E::ScalarField::one()];
        cur = hat_sq;
        for _ in 0..(q-1) {
            powers_of_sq.push(cur);
            cur *= &hat_sq;
        }

        let nu_powers_of_s = powers_of_s.clone()
        .into_par_iter().map(|x| x * &nu)
        .collect::<Vec<_>>();
        let nu_powers_of_sq = powers_of_sq.clone()
        .into_par_iter().map(|x| x * &nu)
        .collect::<Vec<_>>();
        
        let s_g = g_hat.mul(hat_s);
        let sq_h = h_hat.mul(hat_sq);
        let nu_g = g_hat.mul(nu);
        let nu_h = h_hat.mul(nu);
        let hat_s_inv = hat_s.inverse().unwrap();
        let g_0 = g_hat.mul(hat_s_inv);
        let h_0 = h_hat.clone();
        let u = E::pairing(g_0, h_0);
        let tilde_u = E::pairing(
            g_hat.mul(hat_s_inv *  hat_s_inv), h_hat);
        let vec_g = powers_of_s
        .into_par_iter()
        .map(|x| g_hat.mul(x))
        .collect::<Vec<_>>();
        let vec_h = powers_of_sq.clone()
        .into_par_iter()
        .map(|x| h_hat.mul(x))
        .collect::<Vec<_>>();
        let vec_g_prime = powers_of_sq
        .into_par_iter()
        .map(|x| g_hat.mul(x))
        .collect::<Vec<_>>();
        let vec_nu_g = nu_powers_of_s
        .into_par_iter()
        .map(|x| g_hat.mul(x))
        .collect::<Vec<_>>();
        let vec_nu_h = nu_powers_of_sq
        .into_par_iter()
        .map(|x| h_hat.mul(x))
        .collect::<Vec<_>>();
        
        let pp = UniversalParams {
            q,
            num_vars,
            g_hat,
            h_hat,
            s_g,
            sq_h,
            nu_g,
            nu_h,
            g_0,
            h_0,
            u,
            tilde_u,
            vec_g,
            vec_h,
            vec_g_prime,
            vec_nu_g,
            vec_nu_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }

    /// Outputs a commitment to a matrix
    pub fn commit(
        pp: &UniversalParams<E>,
        mat: &Vec<Vec<i32>>,
        hiding_factor: E::ScalarField,
        k: usize,
    ) -> Result<PairingOutput<E>, Error> {

        let m = mat.len();
        let n = mat[0].len();

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial with {:?} coefficients:",
            m * n
        ));

        let vec_g = pp.vec_g[0..m].to_vec();
        let vec_h = pp.vec_h[0..n].to_vec();

        let timer = Instant::now();
        let prepare_g =
            prepare_base_short_g1::<E>(&vec_g, k);
        let mut tier_one_vec = Vec::new();
        for i in 0..m {
            let row = &mat[i];
            let commit_row = msm_g1_short_i32::<E>(
                &prepare_g, &row, k);
            tier_one_vec.push(commit_row);
        }
        let tier1_time = timer.elapsed().as_secs_f64();
        println!("Tier 1 time: {:?}s", tier1_time);

        let result = inner_pairing_product(
            &tier_one_vec, &vec_h)
            + pp.tilde_u.mul(hiding_factor);

        end_timer!(commit_time);
        Ok(result)
    }

    // /// Compute witness polynomial.
    // ///
    // /// The witness polynomial $w(x)$ the quotient of the division (p(x) - p(z)) / (x - z)
    // /// Observe that this quotient does not change with $z$ because
    // /// $p(z)$ is the remainder term. We can therefore omit $p(z)$ when computing the quotient.
    // pub fn compute_witness_polynomial(
    //     p: &P,
    //     point: P::Point,
    //     randomness: &Randomness<E::ScalarField, P>,
    // ) -> Result<(P, Option<P>), Error> {
    //     let divisor = P::from_coefficients_vec(vec![-point, E::ScalarField::one()]);

    //     let witness_time = start_timer!(|| "Computing witness polynomial");
    //     let witness_polynomial = p / &divisor;
    //     end_timer!(witness_time);

    //     let random_witness_polynomial = if randomness.is_hiding() {
    //         let random_p = &randomness.blinding_polynomial;

    //         let witness_time = start_timer!(|| "Computing random witness polynomial");
    //         let random_witness_polynomial = random_p / &divisor;
    //         end_timer!(witness_time);
    //         Some(random_witness_polynomial)
    //     } else {
    //         None
    //     };

    //     Ok((witness_polynomial, random_witness_polynomial))
    // }

    // /// Yields a [`Proof`] with a witness polynomial.
    // pub fn open_with_witness_polynomial<'a>(
    //     powers: &Powers<E>,
    //     point: P::Point,
    //     randomness: &Randomness<E::ScalarField, P>,
    //     witness_polynomial: &P,
    //     hiding_witness_polynomial: Option<&P>,
    // ) -> Result<Proof<E>, Error> {
    //     Self::check_degree_is_too_large(witness_polynomial.degree(), powers.size())?;
    //     let (num_leading_zeros, witness_coeffs) =
    //         skip_leading_zeros_and_convert_to_bigints(witness_polynomial);

    //     let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomial");
    //     let mut w = <E::G1 as VariableBaseMSM>::msm_bigint(
    //         &powers.powers_of_g[num_leading_zeros..],
    //         &witness_coeffs,
    //     );
    //     end_timer!(witness_comm_time);

    //     let random_v = if let Some(hiding_witness_polynomial) = hiding_witness_polynomial {
    //         let blinding_p = &randomness.blinding_polynomial;
    //         let blinding_eval_time = start_timer!(|| "Evaluating random polynomial");
    //         let blinding_evaluation = blinding_p.evaluate(&point);
    //         end_timer!(blinding_eval_time);

    //         let random_witness_coeffs = convert_to_bigints(&hiding_witness_polynomial.coeffs());
    //         let witness_comm_time =
    //             start_timer!(|| "Computing commitment to random witness polynomial");
    //         w += &<E::G1 as VariableBaseMSM>::msm_bigint(
    //             &powers.powers_of_gamma_g,
    //             &random_witness_coeffs,
    //         );
    //         end_timer!(witness_comm_time);
    //         Some(blinding_evaluation)
    //     } else {
    //         None
    //     };

    //     Ok(Proof {
    //         w: w.into_affine(),
    //         random_v,
    //     })
    // }

    // /// On input a polynomial `p` and a `point`, outputs a [`Proof`] for the same.
    // pub fn open<'a>(
    //     powers: &Powers<E>,
    //     p: &P,
    //     point: P::Point,
    //     rand: &Randomness<E::ScalarField, P>,
    // ) -> Result<Proof<E>, Error> {
    //     Self::check_degree_is_too_large(p.degree(), powers.size())?;
    //     let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));

    //     let witness_time = start_timer!(|| "Computing witness polynomials");
    //     let (witness_poly, hiding_witness_poly) = Self::compute_witness_polynomial(p, point, rand)?;
    //     end_timer!(witness_time);

    //     let proof = Self::open_with_witness_polynomial(
    //         powers,
    //         point,
    //         rand,
    //         &witness_poly,
    //         hiding_witness_poly.as_ref(),
    //     );

    //     end_timer!(open_time);
    //     proof
    // }

    // /// Verifies that `value` is the evaluation at `point` of the polynomial
    // /// committed inside `comm`.
    // pub fn check(
    //     vk: &VerifierKey<E>,
    //     comm: &Commitment<E>,
    //     point: E::ScalarField,
    //     value: E::ScalarField,
    //     proof: &Proof<E>,
    // ) -> Result<bool, Error> {
    //     let check_time = start_timer!(|| "Checking evaluation");
    //     let mut inner = comm.0.into_group() - &vk.g.mul(value);
    //     if let Some(random_v) = proof.random_v {
    //         inner -= &vk.gamma_g.mul(random_v);
    //     }
    //     let lhs = E::pairing(inner, vk.h);

    //     let inner = vk.beta_h.into_group() - &vk.h.mul(point);
    //     let rhs = E::pairing(proof.w, inner);

    //     end_timer!(check_time, || format!("Result: {}", lhs == rhs));
    //     Ok(lhs == rhs)
    // }

    // /// Check that each `proof_i` in `proofs` is a valid proof of evaluation for
    // /// `commitment_i` at `point_i`.
    // pub fn batch_check<R: RngCore>(
    //     vk: &VerifierKey<E>,
    //     commitments: &[Commitment<E>],
    //     points: &[E::ScalarField],
    //     values: &[E::ScalarField],
    //     proofs: &[Proof<E>],
    //     rng: &mut R,
    // ) -> Result<bool, Error> {
    //     let check_time =
    //         start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));

    //     let mut total_c = <E::G1>::zero();
    //     let mut total_w = <E::G1>::zero();

    //     let combination_time = start_timer!(|| "Combining commitments and proofs");
    //     let mut randomizer = E::ScalarField::one();
    //     // Instead of multiplying g and gamma_g in each turn, we simply accumulate
    //     // their coefficients and perform a final multiplication at the end.
    //     let mut g_multiplier = E::ScalarField::zero();
    //     let mut gamma_g_multiplier = E::ScalarField::zero();
    //     for (((c, z), v), proof) in commitments.iter().zip(points).zip(values).zip(proofs) {
    //         let w = proof.w;
    //         let mut temp = w.mul(*z);
    //         temp += &c.0;
    //         let c = temp;
    //         g_multiplier += &(randomizer * v);
    //         if let Some(random_v) = proof.random_v {
    //             gamma_g_multiplier += &(randomizer * &random_v);
    //         }
    //         total_c += &c.mul(randomizer);
    //         total_w += &w.mul(randomizer);
    //         // We don't need to sample randomizers from the full field,
    //         // only from 128-bit strings.
    //         randomizer = u128::rand(rng).into();
    //     }
    //     total_c -= &vk.g.mul(g_multiplier);
    //     total_c -= &vk.gamma_g.mul(gamma_g_multiplier);
    //     end_timer!(combination_time);

    //     let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
    //     let affine_points = E::G1::normalize_batch(&[-total_w, total_c]);
    //     let (total_w, total_c) = (affine_points[0], affine_points[1]);
    //     end_timer!(to_affine_time);

    //     let pairing_time = start_timer!(|| "Performing product of pairings");
    //     let result = E::multi_pairing(
    //         [total_w, total_c],
    //         [vk.prepared_beta_h.clone(), vk.prepared_h.clone()],
    //     )
    //     .0
    //     .is_one();
    //     end_timer!(pairing_time);
    //     end_timer!(check_time, || format!("Result: {}", result));
    //     Ok(result)
    // }

    // pub(crate) fn check_degree_is_too_large(degree: usize, num_powers: usize) -> Result<(), Error> {
    //     let num_coefficients = degree + 1;
    //     if num_coefficients > num_powers {
    //         Err(Error::TooManyCoefficients {
    //             num_coefficients,
    //             num_powers,
    //         })
    //     } else {
    //         Ok(())
    //     }
    // }

    // pub(crate) fn check_hiding_bound(
    //     hiding_poly_degree: usize,
    //     num_powers: usize,
    // ) -> Result<(), Error> {
    //     if hiding_poly_degree == 0 {
    //         Err(Error::HidingBoundIsZero)
    //     } else if hiding_poly_degree >= num_powers {
    //         // The above check uses `>=` because committing to a hiding poly with
    //         // degree `hiding_poly_degree` requires `hiding_poly_degree + 1`
    //         // powers.
    //         Err(Error::HidingBoundToolarge {
    //             hiding_poly_degree,
    //             num_powers,
    //         })
    //     } else {
    //         Ok(())
    //     }
    // }

    // pub(crate) fn check_degrees_and_bounds<'a>(
    //     supported_degree: usize,
    //     max_degree: usize,
    //     enforced_degree_bounds: Option<&[usize]>,
    //     p: &'a LabeledPolynomial<E::ScalarField, P>,
    // ) -> Result<(), Error> {
    //     if let Some(bound) = p.degree_bound() {
    //         let enforced_degree_bounds =
    //             enforced_degree_bounds.ok_or(Error::UnsupportedDegreeBound(bound))?;

    //         if enforced_degree_bounds.binary_search(&bound).is_err() {
    //             Err(Error::UnsupportedDegreeBound(bound))
    //         } else if bound < p.degree() || bound > max_degree {
    //             return Err(Error::IncorrectDegreeBound {
    //                 poly_degree: p.degree(),
    //                 degree_bound: p.degree_bound().unwrap(),
    //                 supported_degree,
    //                 label: p.label().to_string(),
    //             });
    //         } else {
    //             Ok(())
    //         }
    //     } else {
    //         Ok(())
    //     }
    // }
}


/// Experiment for SMART PC
pub fn experiment_smart_template<E>(num_vars: usize)
where
    E: Pairing,
{

    let rng = &mut ark_std::test_rng();

    let start_setup = Instant::now();
    let pp = SmartPC::<E>::setup(num_vars, rng).unwrap();
    let setup_time = start_setup.elapsed().as_secs_f64();
    println!("Setup time: {:?}s", setup_time);

    let n : usize = 1 << num_vars/2;
    let k = 8;
    let mat = (0..n).into_par_iter()
    .map(|_|{
        use rand::Rng;
        let mut rng = thread_rng();
        (0..n).into_iter().map(|_|{
            rng.gen_range(-127..127)
        }).collect()
    }).collect();

    let hiding_factor = E::ScalarField::rand(rng);

    let start_commit = Instant::now();
    let comm =
        SmartPC::<E>::commit(&pp, &mat, hiding_factor, k).unwrap();
    let commit_time = start_commit.elapsed().as_secs_f64();
    println!("Commit time: {:?}s", commit_time);

    let mut commit_writer = Vec::new();
    comm.serialize_compressed(&mut commit_writer).unwrap();
    let commit_size = commit_writer.len();
    println!("Commit size: {:?}B", commit_size);


    // let mut commit_writer = Vec::new();
    // comm.serialize_compressed(&mut commit_writer).unwrap();
    // let commit_size = commit_writer.len();
    // println!("Commit size: {:?}B", commit_size);

    // let start_open = Instant::now();
    // let point = E::ScalarField::rand(rng);
    // let value = p.evaluate(&point);
    // let proof = KZG10::<E, P>::open(&ck, &p, point, &rand).unwrap();
    // let open_time = start_open.elapsed().as_secs_f64();
    // println!("Open time: {:?}s", open_time);

    // let verify_time = Instant::now();
    // let verified = KZG10::<E, P>::check(&vk, &comm, point, value, &proof).unwrap();
    // let verify_time = verify_time.elapsed().as_secs_f64()*1000.0;
    // println!("Verify time: {:?}ms", verify_time);

    // let mut writer = Vec::new();
    // proof.serialize_compressed(&mut writer).unwrap();
    // let proof_size = writer.len();
    // println!("Proof size: {:?}B", proof_size);

    // println!(
    //     "KZG10 proof verified: {}",
    //     KZG10::<E, P>::check(&vk, &comm, point, value, &proof).unwrap()
    // );
}

/// Test existing utils function
pub fn test_smart<E:Pairing>() {
    test_utils::<E>();
    // let mut rng = ark_std::test_rng();
    // let degree = 10;
    // let num_vars = 5;
    // let num_samples = 1;
    // let num_queries = 1;
    // let num_coefficients
}