use std::ops::Index;

use crate::{grain::Grain, matrix::Matrix};
use halo2curves::FieldExt;

/// `State` is structure `T` sized field elements that are subjected to
/// permutation
#[derive(Clone, Debug, PartialEq)]
pub struct State<F: FieldExt, const T: usize>(pub(crate) [F; T]);

impl<F: FieldExt, const T: usize> Default for State<F, T> {
    /// The capacity value is 2**64 + (o − 1) where o the output length.
    fn default() -> Self {
        let mut state = [F::zero(); T];
        state[0] = F::from_u128(1 << 64);
        State(state)
    }
}

impl<F: FieldExt, const T: usize> State<F, T> {
    /// Applies sbox for all elements of the state.
    /// Only supports `alpha = 5` sbox case.
    pub(crate) fn sbox_full(&mut self) {
        for e in self.0.iter_mut() {
            let tmp = e.mul(*e);
            e.mul_assign(tmp);
            e.mul_assign(tmp);
        }
    }

    /// Partial round sbox applies sbox to the first element of the state.
    /// Only supports `alpha = 5` sbox case
    pub(crate) fn sbox_part(&mut self) {
        let tmp = self.0[0].mul(self.0[0]);
        self.0[0].mul_assign(tmp);
        self.0[0].mul_assign(tmp);
    }

    /// Adds constants to all elements of the state
    pub(crate) fn add_constants(&mut self, constants: &[F; T]) {
        for (e, constant) in self.0.iter_mut().zip(constants.iter()) {
            e.add_assign(constant)
        }
    }

    /// Only adds a constant to the first element of the state.It is used with
    /// optimized rounds constants where only single element is added in
    /// each partial round
    pub(crate) fn add_constant(&mut self, constant: &F) {
        self.0[0].add_assign(constant)
    }

    /// Copies elements of the state
    pub fn words(&self) -> [F; T] {
        self.0
    }

    /// Second element of the state is the result
    pub(crate) fn result(&self) -> F {
        self.0[1]
    }
}

/// `Spec` holds construction parameters as well as constants that are used in
/// permutation step. Constants are planned to be hardcoded once transcript
/// design matures. Number of partial rounds can be deriven from number of
/// constants.
#[derive(Debug, Clone)]
pub struct Spec<F: FieldExt, const T: usize, const RATE: usize> {
    pub(crate) r_f: usize,
    pub(crate) mds_matrices: MDSMatrices<F, T, RATE>,
    pub(crate) constants: OptimizedConstants<F, T>,
}

impl<F: FieldExt, const T: usize, const RATE: usize> Spec<F, T, RATE> {
    /// Number of full rounds
    pub fn r_f(&self) -> usize {
        self.r_f
    }
    /// Set of MDS Matrices used in permutation line
    pub fn mds_matrices(&self) -> &MDSMatrices<F, T, RATE> {
        &self.mds_matrices
    }
    /// Optimised round constants
    pub fn constants(&self) -> &OptimizedConstants<F, T> {
        &self.constants
    }
}

/// `OptimizedConstants` has round constants that are added each round. While
/// full rounds has T sized constants there is a single constant for each
/// partial round
#[derive(Debug, Clone)]
pub struct OptimizedConstants<F: FieldExt, const T: usize> {
    pub(crate) start: Vec<[F; T]>,
    pub(crate) partial: Vec<F>,
    pub(crate) end: Vec<[F; T]>,
}

impl<F: FieldExt, const T: usize> OptimizedConstants<F, T> {
    /// Returns rounds constants for first part of full rounds
    pub fn start(&self) -> &Vec<[F; T]> {
        &self.start
    }

    /// Returns rounds constants for partial rounds
    pub fn partial(&self) -> &Vec<F> {
        &self.partial
    }

    /// Returns rounds constants for second part of full rounds
    pub fn end(&self) -> &Vec<[F; T]> {
        &self.end
    }
}

/// `MDSMatrices` holds the MDS matrix for the external/internal rounds
#[derive(Debug, Clone)]
pub struct MDSMatrices<F: FieldExt, const T: usize, const RATE: usize> {
    pub(crate) mds_external: MDSMatrix<F, T, RATE>,
    pub(crate) mds_internal: MDSMatrix<F, T, RATE>,
}

impl<F: FieldExt, const T: usize, const RATE: usize> MDSMatrices<F, T, RATE> {
    /// Returns the external MDS matrix
    pub fn mds_external(&self) -> &MDSMatrix<F, T, RATE> {
        &self.mds_external
    }

    /// Returns the internal MDS matrix
    pub fn mds_internal(&self) -> &MDSMatrix<F, T, RATE> {
        &self.mds_internal
    }
}

/// `MDSMatrix` is applied to `State` to achive linear layer of Poseidon
#[derive(Clone, Debug)]
pub struct MDSMatrix<F: FieldExt, const T: usize, const RATE: usize>(pub(crate) Matrix<F, T>);

impl<F: FieldExt, const T: usize, const RATE: usize> Index<usize> for MDSMatrix<F, T, RATE> {
    type Output = [F; T];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.0.0[idx]
    }
}

impl<F: FieldExt, const T: usize, const RATE: usize> MDSMatrix<F, T, RATE> {
    /// Applies `MDSMatrix` to the state
    pub(crate) fn apply(&self, state: &mut State<F, T>) {
        state.0 = self.0.mul_vector(&state.0);
    }

    /// Returns rows of the MDS matrix
    pub fn rows(&self) -> [[F; T]; T] {
        self.0 .0
    }
}


impl<F: FieldExt, const T: usize, const RATE: usize> Spec<F, T, RATE> {
    /// Given number of round parameters constructs new Posedion instance
    /// calculating unoptimized round constants with reference `Grain` then
    /// calculates optimized constants and sparse matrices
    pub fn new(r_f: usize, r_p: usize) -> Self {
        assert_eq!(r_f % 2, 0);
        let (unoptimized_constants, mds_external, mds_internal) = Grain::generate(r_f, r_p);
        let constants = Self::calculate_optimized_constants(r_f, r_p, unoptimized_constants);

        Self {
            r_f,
            constants,
            mds_matrices: MDSMatrices {
                mds_external,
                mds_internal
            },
        }
    }

    // doesn't do any optimization because it's default poseidon2
    fn calculate_optimized_constants(
        r_f: usize,
        r_p: usize,
        constants: Vec<[F; T]>,
    ) -> OptimizedConstants<F, T> {
        let (number_of_rounds, r_f_half) = (r_f + r_p, r_f / 2);
        assert_eq!(constants.len(), number_of_rounds);

        // Calculate optimized constants for first half of the full rounds
        let constants_start: Vec<[F; T]> = constants[0..r_f/2].to_vec();
        let constants_end: Vec<[F; T]> = constants[r_f/2+r_p..].to_vec();

        let mut constants_partial = vec![F::zero(); r_p];
        for i in 0..r_p {
            constants_partial[i] = constants[i + r_f_half][0];
            for j in 1..T {
                assert_eq!(constants[i + r_f_half][j], F::from(0));
            }
        }

        OptimizedConstants {
            start: constants_start,
            partial: constants_partial,
            end: constants_end,
        }
    }
}

#[cfg(test)]
pub(super) mod tests {
    use halo2curves::FieldExt;

    use super::MDSMatrix;
    use crate::grain::Grain;

    /// We want to keep unoptimized parameters to cross test with optimized one
    pub(crate) struct SpecRef<F: FieldExt, const T: usize, const RATE: usize> {
        pub(crate) r_f: usize,
        pub(crate) r_p: usize,
        pub(crate) mds_external: MDSMatrix<F, T, RATE>,
        pub(crate) mds_internal: MDSMatrix<F, T, RATE>,
        pub(crate) constants: Vec<[F; T]>,
    }

    impl<F: FieldExt, const T: usize, const RATE: usize> SpecRef<F, T, RATE> {
        pub(crate) fn new(r_f: usize, r_p: usize) -> Self {
            let (constants, mds_external, mds_internal) = Grain::generate(r_f, r_p);

            /*
            for constant in &constants {
                println!("{:?}", constant);
            }
            */

            SpecRef {
                r_f,
                r_p,
                mds_external,
                mds_internal,
                constants,
            }
        }
    }
}
