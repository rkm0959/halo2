use crate::{Spec, State};
use halo2curves::FieldExt;

/// Poseidon hasher that maintains state and inputs and yields single element
/// output when desired
#[derive(Debug, Clone)]
pub struct Poseidon2<F: FieldExt, const T: usize, const RATE: usize> {
    state: State<F, T>,
    spec: Spec<F, T, RATE>,
    absorbing: Vec<F>,
}

impl<F: FieldExt, const T: usize, const RATE: usize> Poseidon2<F, T, RATE> {
    /// Constructs a clear state poseidon instance
    pub fn new(r_f: usize, r_p: usize) -> Self {
        Self {
            spec: Spec::new(r_f, r_p),
            state: State::default(),
            absorbing: Vec::new(),
        }
    }

    /// Appends elements to the absorption line updates state while `RATE` is
    /// full
    pub fn update(&mut self, elements: &[F]) {
        let mut input_elements = self.absorbing.clone();
        input_elements.extend_from_slice(elements);

        for chunk in input_elements.chunks(RATE) {
            if chunk.len() < RATE {
                // Must be the last iteration of this update. Feed unpermutaed inputs to the
                // absorbation line
                self.absorbing = chunk.to_vec();
            } else {
                // Add new chunk of inputs for the next permutation cycle.
                for (input_element, state) in chunk.iter().zip(self.state.0.iter_mut().skip(1)) {
                    state.add_assign(input_element);
                }
                // Perform intermediate permutation
                self.spec.permute(&mut self.state);
                // Flush the absorption line
                self.absorbing.clear();
            }
        }
    }

    /// Results a single element by absorbing already added inputs
    pub fn squeeze(&mut self) -> F {
        let mut last_chunk = self.absorbing.clone();
        {
            // Expect padding offset to be in [0, RATE)
            debug_assert!(last_chunk.len() < RATE);
        }
        // Add the finishing sign of the variable length hashing. Note that this mut
        // also apply when absorbing line is empty
        last_chunk.push(F::one());
        // Add the last chunk of inputs to the state for the final permutation cycle

        for (input_element, state) in last_chunk.iter().zip(self.state.0.iter_mut().skip(1)) {
            state.add_assign(input_element);
        }

        // Perform final permutation
        self.spec.permute(&mut self.state);
        // Flush the absorption line
        self.absorbing.clear();
        // Returns the challenge while preserving internal state
        self.state.result()
    }
}
