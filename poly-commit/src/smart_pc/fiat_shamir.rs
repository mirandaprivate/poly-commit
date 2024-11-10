//! Utility functions for the adaptive Fiat-Shamir transformation
//! 
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use sha2::{Sha256, Digest};

/// Generate random field element
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
)]
pub struct FiatShamir {
    /// The hash function used by the Fiat-Shamir transform
    pub hasher: Sha256,
    /// The current state of hash value
    pub state: [u8; 32],
}

impl FiatShamir {
    /// Create a new Fiat-Shamir transform
    pub fn new() -> Self {
        let hasher = Sha256::new();
        let state = [0u8; 32];
        Self { hasher, state }
    }

    /// Update the state of the Fiat-Shamir transform
    pub fn push <T>(&mut self, data: &T)
    where
        T: CanonicalSerialize,
    {
        let mut data_bytes = Vec::new();
        data.serialize_compressed(&mut data_bytes).unwrap();
        self.hasher.update(data_bytes);
        self.state.copy_from_slice(self.hasher.finalize_reset().as_slice());
    }

    /// Get the current state of the Fiat-Shamir transform
    pub fn gen_challenge <E:Pairing> (&mut self) -> <E as Pairing>::ScalarField {
        let challenge_bytes = &self.state;
        let challenge =
        E::ScalarField::from_be_bytes_mod_order(challenge_bytes.as_slice());

        self.push::<E::ScalarField>(&challenge);

        challenge
    }

    /// Get the current state of the Fiat-Shamir transform
    pub fn get_state(&self) -> [u8; 32] {
        self.state
    }

}