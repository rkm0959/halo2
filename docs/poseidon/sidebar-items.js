window.SIDEBAR_ITEMS = {"struct":[["MDSMatrices","`MDSMatrices` holds the MDS matrix as well as transition matrix which is also called `pre_sparse_mds` and sparse matrices that enables us to reduce number of multiplications in apply MDS step"],["MDSMatrix","`MDSMatrix` is applied to `State` to achive linear layer of Poseidon"],["Poseidon","Poseidon hasher that maintains state and inputs and yields single element output when desired"],["SparseMDSMatrix","`SparseMDSMatrix` are in `[row], [hat | identity]` form and used in linear layer of partial rounds instead of the original MDS"],["Spec","`Spec` holds construction parameters as well as constants that are used in permutation step. Constants are planned to be hardcoded once transcript design matures. Number of partial rounds can be deriven from number of constants."],["State","`State` is structure `T` sized field elements that are subjected to permutation"]]};