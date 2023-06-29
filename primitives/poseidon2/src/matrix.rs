//! Most of these operations here are not suitable for general purpose matrix
//! operations. Besides vector multiplication other operations are presented
//! with the intention of construction of parameters and are not used in the
//! actual permutation process.

use halo2curves::FieldExt;

#[derive(PartialEq, Debug, Clone)]
pub(crate) struct Matrix<F: FieldExt, const T: usize>(pub(crate) [[F; T]; T]);

impl<F: FieldExt, const T: usize> Default for Matrix<F, T> {
    fn default() -> Self {
        Matrix([[F::zero(); T]; T])
    }
}

impl<F: FieldExt, const T: usize> Matrix<F, T> {
    pub(crate) fn mul_vector(&self, v: &[F; T]) -> [F; T] {
        let mut result = [F::zero(); T];
        for (row, cell) in self.0.iter().zip(result.iter_mut()) {
            for (a_i, v_i) in row.iter().zip(v.iter()) {
                *cell += *v_i * *a_i;
            }
        }
        result
    }
}
