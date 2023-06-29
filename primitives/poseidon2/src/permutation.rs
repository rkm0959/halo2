use crate::spec::{Spec, State};
use halo2curves::FieldExt;

impl<F: FieldExt, const T: usize, const RATE: usize> Spec<F, T, RATE> {
    /// Applies the Poseidon2 permutation to the given state
    pub fn permute(&self, state: &mut State<F, T>) {
        let (mds_external, mds_internal) = (self.mds_matrices().mds_external(), self.mds_matrices().mds_internal());

        mds_external.apply(state);

        for constants in self.constants.start().iter() {
            state.add_constants(constants);
            state.sbox_full();
            mds_external.apply(state);
        }

        for constant in self.constants.partial().iter() {
            state.add_constant(constant); 
            state.sbox_part();
            mds_internal.apply(state);
        }

        for constants in self.constants.end().iter() {
            state.add_constants(constants);
            state.sbox_full();
            mds_external.apply(state);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::State;
    use crate::spec::{tests::SpecRef, Spec};
    use group::ff::PrimeField;
    use halo2curves::bn256::Fr;
    use halo2curves::FieldExt;

    /// We want to keep unoptimized poseidion construction and permutation to
    /// cross test with optimized one
    impl<F: FieldExt, const T: usize, const RATE: usize> SpecRef<F, T, RATE> {
        fn permute(&self, state: &mut State<F, T>) {
            let (r_f, r_p) = (self.r_f / 2, self.r_p);

            self.mds_external.apply(state);

            for constants in self.constants.iter().take(r_f) {
                state.add_constants(constants);
                state.sbox_full();
                self.mds_external.apply(state);
            }

            for constants in self.constants.iter().skip(r_f).take(r_p) {
                state.add_constants(constants); // only one here, but no need to optimize here
                state.sbox_part();
                self.mds_internal.apply(state);
            }

            for constants in self.constants.iter().skip(r_f + r_p) {
                state.add_constants(constants);
                state.sbox_full();
                self.mds_external.apply(state);
            }
        }
    }

    #[test]
    fn cross_test() {
        use halo2curves::group::ff::Field;
        use rand_core::OsRng;
        use std::time::Instant;

        macro_rules! run_test {
            (
                $([$RF:expr, $RP:expr, $T:expr, $RATE:expr]),*
            ) => {
                $(
                    {
                        const R_F: usize = $RF;
                        const R_P: usize = $RP;
                        const T: usize = $T;
                        const RATE: usize = $RATE;
                        let mut state = State(
                            (0..T)
                                .map(|_| Fr::random(OsRng))
                                .collect::<Vec<Fr>>()
                                .try_into().unwrap(),
                        );
                        let spec = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
                        let mut state_expected = state.clone();
                        spec.permute(&mut state_expected);

                        let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
                        let now = Instant::now();
                        {
                            spec.permute(&mut state);
                        }
                        let elapsed = now.elapsed();
                        println!("Elapsed: {:.2?}", elapsed);
                        assert_eq!(state_expected, state);
                    }
                )*
            };
        }
        run_test!([8, 56, 3, 2]);
        run_test!([8, 56, 4, 3]);
        run_test!([8, 57, 8, 7]);
        run_test!([8, 57, 12, 11]);
    }

    #[test]
    fn test_against_test_vectors() {
        // from the vectors in https://github.com/HorizenLabs/poseidon2/blob/main/poseidon2_rust_params.sage
        {
            const R_F: usize = 8;
            const R_P: usize = 56;
            const T: usize = 3;
            const RATE: usize = 2;

            let state = State(
                vec![0u64, 1, 2]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "21882471761025344482456282050943515707267606647948403374880378562101343146243",
                "9030699330013392132529464674294378792132780497765201297316864012141442630280",
                "9137931384593657624554037900714196568304064431583163402259937475584578975855",
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
        {
            const R_F: usize = 8;
            const R_P: usize = 56;
            const T: usize = 3;
            const RATE: usize = 2;

            let state = State(
                vec![Fr::from(2).pow(&[64, 0, 0, 0]), Fr::from(1), Fr::from(2)]
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "21658538694766004233136855742940871205577528902229351089041313470429025882597",
                "19296375964009402039404592001274736735894410486509328719258134543103656624902",
                "11640765419168198000719743597189611575188281354497517721813259831478696667706",
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
        {
            const R_F: usize = 8;
            const R_P: usize = 56;
            const T: usize = 4;
            const RATE: usize = 3;

            let state = State(
                vec![0u64, 1, 2, 3]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "20662024781356425604712282113410366519292003656984552908876205683335846693377",
                "1038338918110681223973284474158620201510921073920677497393589062792853859568",
                "17247817174748483971912051526215045751313667531175946679417656805949142294500",
                "13811986223557379116635132669359609086888098128791671082914978079180962867260",
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
        {
            const R_F: usize = 8;
            const R_P: usize = 57;
            const T: usize = 8;
            const RATE: usize = 7;

            let state = State(
                vec![0u64, 1, 2, 3, 4, 5, 6, 7]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "17996387074513178554865200666798625375796202432664131140470595843161387418754", "6519623664793092174849445745730156549937019721325081631911441035361215715490", "15019573654943701764103754286888831501674271555506365891183749798590184995997", "18922559970423556479308188907943183991206778671083370007884963708379149903808", "18788579285582264519982507685282888088815173608105137299326524527726935412283", "9910189047954102661591142129580177705649591816956040874337745316431043262850", "19015531005769431747078148486509559370200972647172060755681977331792915725441", "19553857394421899094377801204336606079539593824011190383971537103654264732633"
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
        {
            const R_F: usize = 8;
            const R_P: usize = 57;
            const T: usize = 12;
            const RATE: usize = 11;

            let state = State(
                vec![0u64, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "1745989064378017317504544419661115956248317509917596575176363569643799142334", "6496133573199966431066054171498945278755874003697713663279764444643329651660", "16972936058828011533307967308063106589693202218195273462056534317946884181842", "11928434893972397031467362379699105906495772180101221872275448576925487477436", "11602401632385076528491334043916760602913028641920062889417726651201044434177", "3384137084037436360403166331755652231301809506355712752229800702508814529628", "18574321955856137637299911851191266356408252699912362648355664428714715739270", "5768738003881838537577414013313012972716892681764978615456623304699593400632", "13616584935433632834317064726010376032813675422184858635946334852412802257769", "7224127608005614511560477486685559795094709765884793782106711176158970279563", "15882960620698794140850040250591091226099106674344621794321485459750023397901", "19167632659951096956673065333410962999212009098364495429373842266953566010993"
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }
            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
    }
}
