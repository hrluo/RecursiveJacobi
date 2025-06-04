# Readme File

# Code Organization and Usage
The code has been tested with MATLAB R2024a, on an AMD Ryzen 9 7945HX CPU with Radeon Graphics.

The codes for experiments are organized into separate folders, each named according to the corresponding figures and tables in the paper. Each folder can be run independently. For experiments that are time-consuming, the precomputed data are provided in corresponding `.mat` files. To reproduce the figures, simply load the `.mat` file and run the plotting sections in the provided scripts.

---
# `Classical Jacobi`

## Description
`ClassicalJacobi` is a flexible implementation of the Classical Jacobi method that supports multiple strategies for computing rotation angles.

## Inputs
- **`A`**: Symmetric matrix to be diagonalized.
- **`eps_threshold`**: Tolerance for convergence, based on off-diagonal norms.
- **`method`**: Specifies the rotation strategy. Options are `'eig'`, `'trig'` and `'adversarial'`:
  - `'eig'`: Uses MATLAB's buit-in `eig()` function on the 2x2 submatrices.
  - `'trig'`: Computes rotation angles explicitly using trigonometric formulas.
  - `'adversarial'`: Deliberately adds pi/2 to the rotation angle computed by `trig` to hinder convergence, following [*The Cyclic Jacobi Method for Computing the Principal Values of a Complex Matrix*](https://doi.org/10.1090/S0002-9947-1960-0109825-2), by Forsythe and Henrici.
- **`n3_ratio`**: Maximum computational budget as a fraction of `n^3` floating-point operations.

## Outputs
- **`Q`**: Orthogonal matrix of eigenvectors.
- **`D`**: Diagonal matrix of eigenvalues.
- **`flops`**: Total floating-point operations performed.
- **`sweeps`**: Total number of sweeps performed
- **`sweep_OffNorm_history`**: History recorded at the end of each sweep in the format `[flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm]`.

## Algorithm Details
1. **Initialization**:
   - Validate matrix symmetry.
   - Initialize `Q` as the identity matrix and set `flops`, `sweeps` to zero.
   - Start tracking execution history as `sweep_OffNorm_history`.

2. **Iterative Updates**:
   - Depending on the `method`:
     - **`'eig'`**: Diagonalize 2x2 submatrices using eigen decomposition.
     - **`'trig'`**: Use trigonometric formulas to compute rotation angles.
     - **`'adversarial'`**: Adding pi/2 to the trigonometric rotation angle to resist convergence
   - Update matrix `A` and orthogonal matrix `Q` with each rotation.
   - Update cumulative FLOPs and `sweep_OffNorm_history`.

3. **Convergence**:
   - Terminate when the off-diagonal norm falls below `eps_threshold` or the FLOP budget is exhausted.

---

# `Block Jacobi`

## Description
`BlockJacobi` performs a block Jacobi method for symmetric eigenproblems. The method partitions the input matrix into blocks and iteratively applies orthogonal transformations to reduce off-diagonal blocks. It supports multiple orderings, including row-cyclic, column-cyclic, and random ordering, as well as various pivoting strategies, including QR decomposition with column pivoting (QRCP, following [*A Global Convergence Proof for Cyclic Jacobi Methods with Block Rotations*](https://doi.org/10.1137/090748548) by Drmac), LU decomposition with partial pivoting (LUPP), and randomized pivoting.

## Inputs
- **`A`**: Symmetric matrix to be diagonalized.
- **`blockSizes`**: Array specifying the sizes of each block. The sum of `blockSizes` must equal the size of `A`.
- **`eps_threshold`**: Tolerance for convergence, based on off-diagonal norms.
- **`ordering`**: Specifies the iteration order. Options are:
  - `'columncyclic'`: Cyclically iterate over column pairs.
  - `'rowcyclic'`: Cyclically iterate over row pairs.
  - `'rowcyclic'`: Randomly iterate over block pairs.
- **`pivot_method`**: Specifies the pivoting strategy for block transformations. Options are:
  - `'eig'`: Directly uses eigen-decomposition for block transformations without any pivoting (i.e. vanilla block Jacobi).
  - `'qrcp'`: Applies QRCP to permute the rotation matrix in addition to solving the block subproblem, which guarantees convergence of block Jacobi.
  - `'lupp'`: Applies LUPP to permute the rotation matrix in addition to solving the block subproblem, which guarantees convergence of block Jacobi.
  - `'fastlu'`: Replaces MATLAB's built-in `lu()` function with a recursive LU decomposition with partial pivoting, which requires fewer FLOPs for galactic matrices.
  - `'random'`: Applies random permutation to the rotation matrix, which serves as a comparison to QRCP and LUPP
- **`n3_ratio`**: Maximum computational budget as a fraction of `n^3` floating-point operations.

## Outputs
- **`Q`**: Orthogonal matrix of eigenvectors.
- **`D`**: Diagonal matrix of eigenvalues.
- **`flops`**: Total floating-point operations performed.
- **`sweeps`**: Total number of sweeps performed
- **`sweep_OffNorm_history`**: History recorded at the end of each sweep in the format `[flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm]`.

## Algorithm Details
1. **Initialization**:
   - Validate matrix symmetry.
   - Initialize `Q` as the identity matrix and set `flops`, `sweeps` to zero.
   - Start tracking execution history as `sweep_OffNorm_history`.

2. **Iterative Updates**:
   - Loop over block pairs following the specified `ordering` (`'rowcyclic'`, `'columncyclic'`, or `'random'`).
   - For each block pair:
   
      - Compute the rotation matrix based on the specified `pivot_method`:
        - `'eig'`: Use eigenvalue decomposition for the selected block pair.
        - `'qrcp'`: Apply QRCP to permute the rotation matrix after solving the block subproblem.
        - `'lupp'`: Apply LUPP to permute the rotation matrix after solving the block subproblem.
        - `'fastlu'`: Replace MATLAB buiilt-in `lu()` with recursive LU decomposition.
        - `'random'`: Apply random permutation to the rotational matrix.
    
      - Update matrix `A` and orthogonal matrix `Q` after each rotation.
      - Update cumulative FLOPs and `sweep_OffNorm_history`.

3. **Convergence**:
   - Terminate when the off-diagonal norm falls below `eps_threshold` or the FLOP budget is exhausted.

---
# Recursive Jacobi Methods

Below provides detailed descriptions and explanations of the Recursive Jacobi methods implemented in the provided MATLAB files. These methods solve symmetric eigenproblems using various techniques for the subproblem within each block such as LU decomposition, QR decomposition. Each function terminates when the off-diagonal norm falls below `eps_threshold` or the FLOP budget is exhausted.

---

## `RecursiveJacobiplain`

## Description
This function implements the vanilla version of the Recursive Jacobi method for symmetric eigenproblems. It directly applies Algorithm 3 without pivoting (i.e., lines 11 - 13) , which serves as the baseline for all other specialized methods.

## Inputs
- **`A`**: Symmetric matrix to be diagonalized.
- **`n_threshold`**: Threshold for solving the problem as base case directly, based on available fast memory.
- **`f`**: the log block size parameter, strictly between 0 and 1.
- **`eps_threshold`**: Tolerance for convergence, based on off-diagonal norms.
- **`recdepth`**: Current recursion depth
- **`break_flag`**: Flag to halt recursion prematurely.
- **`n3_ratio`**: Maximum computational budget as a fraction of `n^3` floating-point operations.

## Outputs
- **`Q`**: Orthogonal matrix of eigenvectors.
- **`D`**: Diagonal matrix of eigenvalues.
- **`flops`**: Total floating-point operations performed.
- **`sweeps`**: Total number of sweeps performed
- **`sweep_OffNorm_history`**: History recorded at the highest level and at the end of each ``sweep" in the format [flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm].
---
The below variants of Recursive Jacobi methods share similar inputs and outputs parameters compared to the `RecursiveJacobiplain` functions, therefore we only states the main difference:

## Algorithm Details
1. **Initialization**:
   - Validate matrix symmetry.
   - Get input problem size `n` and log block parameter `b = n^f`
   - Initialize `Q` as the identity matrix and set `flops`, `sweeps` to zero.
   - If `rec_depth == 0`, start tracking execution history as `sweep_OffNorm_history`.

2. **Solve bottom case directly**
   - If `n <= n_threshold` or `2b >= n`, solve the problem directly.

3. **Recursive Updates**:
   - Iterate thorugh each sub-problem `A_hat`:
   
      - Solve the sub-problem by recursively calling `RecursiveJacobiplain(A_hat)`
    
      - Update matrix `A`, orthogonal matrix `Q` and cumulative FLOPs after each recursive call.
   - If `rec_depth == 0`, update `sweep_OffNorm_history`.

4. **Convergence**:
   - Terminate when the off-diagonal norm falls below `eps_threshold` or the FLOP budget is exhausted.

## `RecursiveJacobiLUPP`

**Difference**: 
- Implements LUPP in lines 11 - 13 to ensure convergence.

## `RecursiveJacobifastLUPP`
**Difference**: 
- The same as RecursiveJacobiLUPP, but the built-in MATLAB `lu()` function is replaced with recursive LUPP, which achieves the optimal asymptotic complexity for galactic matrices.

## `RecursiveJacobiQRCP`

**Difference**: 
- Implements QRCP in lines 11 - 13 to ensure convergence, however at the expense of higher computational cost.

---
# Jacobi SVD

## Description
`blockOneSidedJacobi` implements the classical one-sided Jacobi method to compute the singular value decomposition (SVD) of a square matrix. The algorithm applies a sequence of orthogonal rotations to jointly diagonalize the right singular vectors, driving the matrix towards orthogonal form.

## Inputs
- **`G`**: n x n real matrix to be factorized (assumed square)
- **`b`**: Block size (must divide `n` exactly)
- **`tol`**: Stopping tolerance for "near orthogonality"
- **`maxSweeps`**: Maximum number of outer loop sweeps
- **`wantVectors`**: Boolean flag indicating whether sigular vectors are needed

## Outputs
- **`U`**: Approximate left singular vectors
- **`Sigma`**: Diagonal matrix of singular values
- **`V`**: Approximate right singular vectors

## Algorithm Details
1. **Initialization**:

   * Check that input matrix `G` is square.
   * Verify that block size `b` divides `n` exactly.
   * Compute number of blocks `nb = n / b`.
   * Initialize `V = eye(n)` if `wantVectors == true`; otherwise set `V = []`.

2. **Outer Loop**:

   * Repeat for at most `maxSweeps` outer sweeps:

     * If `G` is "near orthogonal", terminate early.
     * For each block pair `(I, J)` with `1 <= I < J <= nb`:

       * Extract column indices `colsI` and `colsJ` for blocks `I` and `J`.
       * Form submatrix `Gsub = [G(:, colsI), G(:, colsJ)]`.
       * Compute `A_hat = Gsub' * Gsub`.
       * If `A_hat` is far from digonal:

         * Compute eigen-decomposition `[V_hat, ~] = eig(A_hat)`.
         * Optional: If `det(V_hat) < 0`, flip first column of `V_hat` to ensure determinant positivity.
         * Apply rotation: update `G` accordingly.
         * If `wantVectors == true`, update `V` accordingly.

3. **Postprocessing**:

   * Compute column norms of `G` to form diagonal matrix `Sigma`.
   * If `wantVectors == true`, compute `U = G * inv(Sigma)`; otherwise set `U = []`.

---
# Auxiliary functions

## Convergence Measure

We simultaneously track two types of off-diagonal norm to measure convergence.

### `normOffDiag`
Computes the maximum off-diagonal absolute value of the input matrix.

### `FroNormOffDiag`
Computes the off-diagonal Frobenius norm of the input matrix.

## Adversarial Input Matrices

In addition to standard random generation via `A = randn(n_size); A = (A + A')/2`, three types of adversarial symmetric matrices are provided to evaluate the algorithms under challenging scenarios.

### `generate_nearly_signed_permutation_matrix`
$\mathrm{A} = \mathrm{Q} \mathrm{D} \mathrm{Q}^{T}$ where $\mathrm{D}$ is a standard normal random diagonal matrix, and $\mathrm{Q}$ is obtained from a QR decomposition of $\mathrm{I} + \delta \mathrm{G}$ for $ \mathrm{G} \sim \mathcal{N}(0, 1)^{n \times n}$ and (small) $\delta > 0$, i.e. a small pertubation of the identity matrix.

### `generate_scaled_Hadamard_matrix`
$\mathrm{A} = \mathrm{Q} \mathrm{D} \mathrm{Q}^{T}$ where $\mathrm{D}$ is a standard normal random diagonal matrix, and $\mathrm{Q} = \frac{1}{\sqrt{n}} \cdot \text{Hadamard}(n),n=2^N$, a scaled Hadardmard matrix.

### `generate_spike_spectrum_matrix`
$\mathrm{A} = \mathrm{Q} \mathrm{\Lambda} \mathrm{Q}^{T}$ where $\mathrm{\Lambda}$ is diagonal  and $\mathrm{Q}$ is obtained from the QR decomposition of a random Gaussian matrix. To drive up the condition number of $\mathrm{A}$, the diagonal of $\mathrm{\Lambda}$ contains a pre-specified number of dominant eigenvalues (i.e., spikes) that is at least a ratio larger than other smaller ones.

## `fastLU`
Recursive implementation of LUPP, according to [*Fast Linear Algebra is Stable*](https://doi.org/10.1007/s00211-007-0114-x), by Demmel, Dumitriu and Holtz.


---
## References
**Citation**
If you use our code, please refer to 
[LICENSE](https://github.com/hrluo/RecursiveJacobi/blob/master/LICENSE) and please cite our paper using following BibTeX item (we need to change this once submitted):

    @article{2025recursivejacobi,
        title={Recursive Jacobi Methods},
        author={James W. Demmel, Hengrui Luo, Ryan Schneider, Yifu Wang},
        year={2025},
        eprint={https://arxiv.org/abs/xxxx.yyyy},
        archivePrefix={arXiv},
        primaryClass={math.LA}
    }

Thank you again for the interest and please reach out if you have further questions.
