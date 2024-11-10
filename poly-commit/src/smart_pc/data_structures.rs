/// Defines the data structures used in SMART-PC
// use crate::*;
use ark_ec::pairing::{Pairing,PairingOutput};
    // AdditiveGroup, AffineRepr
// };
// use ark_ff::ToConstraintField;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress,
    SerializationError, Valid,
    Validate,
};
use ark_std::io::{Read, Write};


/// `UniversalParams` are the universal parameters for the SMART-PC scheme.
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct UniversalParams<E: Pairing> {
    /// Maximum supportted matrix dimension
    pub q: usize,
    /// Number of variables
    pub num_vars: usize,
    /// G1 generator
    pub g_hat: E::G1,
    /// G2 generator
    pub h_hat: E::G2,
    /// G1 element
    pub s_h: E::G2,
    /// G2 element
    pub sq_h: E::G2,
    /// G1 element
    pub nu_g: E::G1,
    /// G2 element
    pub nu_h: E::G2,
    /// G1 element
    pub g_0: E::G1,
    /// G2 element
    pub h_0: E::G2,
    /// Gt base
    pub u: PairingOutput<E>,
    /// Gt base for randomness
    pub tilde_u: PairingOutput<E>,
    /// G1 vector for commitment
    pub vec_g: Vec<E::G1>,
    /// G2 vector for commitment
    pub vec_h: Vec<E::G2>,
    /// G1 vector for opening
    pub vec_g_prime: Vec<E::G1>,
    /// G1 vector for opening
    pub vec_nu_g: Vec<E::G1>,
    /// G2 vector for opening
    pub vec_nu_h: Vec<E::G2>,
}



impl<E: Pairing> Valid for UniversalParams<E> {
    fn check(&self) -> Result<(), SerializationError> {
        self.q.check()?;
        self.num_vars.check()?;
        self.g_hat.check()?;
        self.h_hat.check()?;
        self.s_h.check()?;
        self.sq_h.check()?;
        self.nu_g.check()?;
        self.nu_h.check()?;
        self.g_0.check()?;
        self.h_0.check()?;
        self.u.check()?;
        self.tilde_u.check()?;
        self.vec_g.check()?;
        self.vec_h.check()?;
        self.vec_g_prime.check()?;
        self.vec_nu_g.check()?;
        self.vec_nu_h.check()?;
        Ok(())
    }
}


impl<E: Pairing> CanonicalSerialize for UniversalParams<E> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.q.serialize_with_mode(&mut writer, compress)?;
        self.num_vars.serialize_with_mode(&mut writer, compress)?;
        self.g_hat.serialize_with_mode(&mut writer, compress)?;
        self.h_hat.serialize_with_mode(&mut writer, compress)?;
        self.s_h.serialize_with_mode(&mut writer, compress)?;
        self.sq_h.serialize_with_mode(&mut writer, compress)?;
        self.nu_g.serialize_with_mode(&mut writer, compress)?;
        self.nu_h.serialize_with_mode(&mut writer, compress)?;
        self.g_0.serialize_with_mode(&mut writer, compress)?;
        self.h_0.serialize_with_mode(&mut writer, compress)?;
        self.u.serialize_with_mode(&mut writer, compress)?;
        self.tilde_u.serialize_with_mode(&mut writer, compress)?;
        self.vec_g.serialize_with_mode(&mut writer, compress)?;
        self.vec_h.serialize_with_mode(&mut writer, compress)?;
        self.vec_g_prime.serialize_with_mode(&mut writer, compress)?;
        self.vec_nu_g.serialize_with_mode(&mut writer, compress)?;
        self.vec_nu_h.serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.q.serialized_size(compress)
            + self.num_vars.serialized_size(compress)
            + self.g_hat.serialized_size(compress)
            + self.h_hat.serialized_size(compress)
            + self.s_h.serialized_size(compress)
            + self.sq_h.serialized_size(compress)
            + self.nu_g.serialized_size(compress)
            + self.nu_h.serialized_size(compress)
            + self.g_0.serialized_size(compress)
            + self.h_0.serialized_size(compress)
            + self.u.serialized_size(compress)
            + self.tilde_u.serialized_size(compress)
            + self.vec_g.serialized_size(compress)
            + self.vec_h.serialized_size(compress)
            + self.vec_g_prime.serialized_size(compress)
            + self.vec_nu_g.serialized_size(compress)
            + self.vec_nu_h.serialized_size(compress)
    }
}

impl<E: Pairing> CanonicalDeserialize for UniversalParams<E> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let q = usize::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let num_vars = usize::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let g_hat = E::G1::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let h_hat = E::G2::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let s_h = E::G2::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let sq_h = E::G2::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let nu_g = E::G1::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let nu_h = E::G2::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let g_0 = E::G1::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let h_0 = E::G2::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let u = PairingOutput::<E>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let tilde_u = PairingOutput::<E>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let vec_g = Vec::<E::G1>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let vec_h = Vec::<E::G2>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let vec_g_prime = Vec::<E::G1>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let vec_nu_g = Vec::<E::G1>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;
        let vec_nu_h = Vec::<E::G2>::deserialize_with_mode(
            &mut reader, compress, Validate::No)?;

        let result = Self {
            q,
            num_vars,
            g_hat,
            h_hat,
            s_h,
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
        if let Validate::Yes = validate {
            result.check()?;
        }

        Ok(result)
    }
}


/// `Trans` are the transcript for the SMART-PC scheme.
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Trans<E: Pairing> {
    /// The L in Bulletproofs with hiding factors
    pub vec_l_tilde: Vec<PairingOutput<E>>,
    /// The R in Bulletproofs with hiding factors
    pub vec_r_tilde: Vec<PairingOutput<E>>,
    /// Hiding commitment to hat_a
    pub com_rhs_tilde: PairingOutput<E>,
    /// V_G in G1
    pub v_g: E::G1,
    /// V_H in G2
    pub v_h: E::G2,
    /// V_G_prime in G2
    pub v_g_prime: E::G1,
    /// V_H_prime in G1
    pub v_h_prime: E::G2,
    /// W_G in G1
    pub w_g: E::G1,
    /// W_H in G1
    pub w_h: E::G1,
    /// Schnorr 1 Group element
    pub schnorr_1_f: PairingOutput<E>,
    /// Schnorr 1 Zp element
    pub schnorr_1_z: E::ScalarField,
    /// Schnorr 2 Group element
    pub schnorr_2_f: PairingOutput<E>,
    /// Schnorr 2 Zp element 1
    pub schnorr_2_z_1: E::ScalarField,
    /// Schnorr 2 Zp element 2
    pub schnorr_2_z_2: E::ScalarField,
}

impl<E: Pairing> Valid for Trans<E> {
    fn check(&self) -> Result<(), SerializationError> {
        self.vec_l_tilde.check()?;
        self.vec_r_tilde.check()?;
        self.com_rhs_tilde.check()?;
        self.v_g.check()?;
        self.v_h.check()?;
        self.v_g_prime.check()?;
        self.v_h_prime.check()?;
        self.w_g.check()?;
        self.w_h.check()?;
        self.schnorr_1_f.check()?;
        self.schnorr_1_z.check()?;
        self.schnorr_2_f.check()?;
        self.schnorr_2_z_1.check()?;
        self.schnorr_2_z_2.check()?;
        Ok(())
    }
}

impl<E: Pairing> CanonicalSerialize for Trans<E> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.vec_l_tilde.serialize_with_mode(&mut writer, compress)?;
        self.vec_r_tilde.serialize_with_mode(&mut writer, compress)?;
        self.com_rhs_tilde.serialize_with_mode(&mut writer, compress)?;
        self.v_g.serialize_with_mode(&mut writer, compress)?;
        self.v_h.serialize_with_mode(&mut writer, compress)?;
        self.v_g_prime.serialize_with_mode(&mut writer, compress)?;
        self.v_h_prime.serialize_with_mode(&mut writer, compress)?;
        self.w_g.serialize_with_mode(&mut writer, compress)?;
        self.w_h.serialize_with_mode(&mut writer, compress)?;
        self.schnorr_1_f.serialize_with_mode(&mut writer, compress)?;
        self.schnorr_1_z.serialize_with_mode(&mut writer, compress)?;
        self.schnorr_2_f.serialize_with_mode(&mut writer, compress)?;
        self.schnorr_2_z_1.serialize_with_mode(&mut writer, compress)?;
        self.schnorr_2_z_2.serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.vec_l_tilde.serialized_size(compress)
            + self.vec_r_tilde.serialized_size(compress)
            + self.com_rhs_tilde.serialized_size(compress)
            + self.v_g.serialized_size(compress)
            + self.v_h.serialized_size(compress)
            + self.v_g_prime.serialized_size(compress)
            + self.v_h_prime.serialized_size(compress)
            + self.w_g.serialized_size(compress)
            + self.w_h.serialized_size(compress)
            + self.schnorr_1_f.serialized_size(compress)
            + self.schnorr_1_z.serialized_size(compress)
            + self.schnorr_2_f.serialized_size(compress)
            + self.schnorr_2_z_1.serialized_size(compress)
            + self.schnorr_2_z_2.serialized_size(compress)
    }
}