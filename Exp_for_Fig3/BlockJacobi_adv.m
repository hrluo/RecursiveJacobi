% MATLAB Implementation of Block Jacobi Algorithm
function [Q, D, flops, sweeps, sweep_OffNorm_history] = BlockJacobi_adv(A, blockSizes, eps_threshold, ordering, pivot_method, n3_ratio)
    % Block Jacobi Algorithm with adversarial bottom cases for Symmetric Eigenproblem
    % Inputs:
    %   A - symmetric matrix to be diagonalized
    %   blockSizes - sizes of block partitions
    %   eps_threshold - tolerance for off-diagonal norm
    %   ordering - 'columncyclic' 
    %            - 'rowcyclic'
    %            - 'random'
    %   pivot_method - 'eig'          % vanilla block Jacobi
    %                - 'qrcp'         % Based on GLOBAL CONVERGENCE PROOF FOR CYCLIC JACOBI METHODS WITH BLOCK ROTATIONS by ZLATKO DRMAČ, 2009.
    %                - 'lupp'         % block Jacobi with LU with partial pivoting
    %                - 'fastlu'       % recursive implementation of LUPP, based on Fast Linear Algebra is Stable by James Demmel, Ioana Dumitriu, and Olga Holtz, 2007.
    %                - 'random'       % randomly permute the rotational matrix, as a comparison to QRCP and LUPP
    %   n3_ratio - FLOP threshold relative to n^3
    % Outputs:
    %   Q - orthogonal matrix of eigenvectors
    %   D - near-diagonal matrix of eigenvalues
    %   flops - total FLOPs performed
    %   sweeps - total sweeps performed
    %   sweep_OffNorm_history - a log of execution, recorded in the format [flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm], with entries captured at the end of each sweep

    % Validate input
    if ~isequal(A, A')
        error('Input matrix must be Hermitian.');
    end

    if ~any(strcmp(ordering, {'columncyclic', 'rowcyclic', 'random'}))
        error("Pivoting method must be either 'columncyclic', 'rowcyclic', or 'random'.");
    end

    % Initialize variables
    n = size(A, 1);
    m = length(blockSizes);
    indices = mat2cell(1:n, 1, blockSizes);
    Q = eye(n);
    flops = 0;
    sweeps = 0;
    off_diag_norm = normOffDiag(A);
    fro_norm = FroNormOffDiag(A);

    % record history
    sweep_OffNorm_history = [0, 0, off_diag_norm, fro_norm];
    
    % Generate block pairs based on pivot method
    if strcmp(ordering, 'columncyclic')
        blockPairs = generateColumnCyclicPairs(m);
    elseif strcmp(ordering, 'rowcyclic')
        blockPairs = generateRowCyclicPairs(m);
    end

    % Main loop for block Jacobi iterations
    while flops < n3_ratio * n^3
        sweeps = sweeps + 1;
        
        if strcmp(ordering, 'random')
            blockPairs = generateRandomPairs(m);
        end
        
        for pairIdx = 1:size(blockPairs, 1)
            % Select block pair
            i = blockPairs(pairIdx, 1);
            j = blockPairs(pairIdx, 2);

            % Get block indices
            blockIdx_i = indices{i};
            blockIdx_j = indices{j};

            % Extract block submatrix
            rows = [blockIdx_i, blockIdx_j];
            cols = [blockIdx_i, blockIdx_j];
            blockSubMatrix = A(rows, cols);
            blockSubMatrix = (blockSubMatrix + blockSubMatrix') / 2;   % ensure symmetry

            if normOffDiag(blockSubMatrix) >= eps_threshold
                % [U, ~] = eig(blockSubMatrix);
                % flops = flops + (8 + 2/3) * size(blockSubMatrix, 1)^3;
                % Here we substitue the original eig() with self-written adversarial classical Jacobi to resist convergence
                [U, ~, adv_flops, ~, ~] = classicalJacobi(blockSubMatrix, eps_threshold, "adversarial", Inf);
                flops = flops + adv_flops;
                
                % Update the orthogonal matrix based on different methods we choose
                switch lower(pivot_method)
                    case 'eig'
                        % Do nothing. In this case, the algorithm is reduced to the vanilla block Jacobi method.
                    case 'qrcp'
                        % This is the implementation as in the ZLATKO DRMAČ, 2009 paper.
                        b1 = size(blockIdx_i, 2);
                        m1 = size(U, 1);
                        n1 = size(U, 2);
    
                        U_11_12 = U(1:b1, :);
                        [~, ~, P_qrcp] = qr(U_11_12);
                        U = U * P_qrcp;
                        flops_qr = 2 * n1 * b1^2 - (2 / 3) * b1^3;   % Swap n1, b1 since the matrix is fat and short.
                        flops_perm = 2 * m1 * n1^2;
                        flops = flops + flops_qr + flops_perm;
                    case 'lupp'
                        % This is the implementation usnig LUPP instead of QRCP to permute the rotational matrix U
                        b1 = size(blockIdx_i, 2);
                        m1 = size(U, 1);
                        n1 = size(U, 2);
    
                        U_11_12 = U(1:b1, :);
                        [~, ~, P_lupp] = lu(U_11_12');
                        U = U * P_lupp';
                        flops_lu = n1 * b1^2 - (1 / 3) * b1^3;
                        flops_perm = 2 * m1 * n1^2;
                        flops = flops + flops_lu + flops_perm;
                    case 'fastlu'
                        % This is the implementation of LUPP using a recursive algorithm, based on Fast Linear Algebra is Stable by James Demmel, Ioana Dumitriu, and Olga Holtz, 2007.
                        b1 = size(blockIdx_i, 2);
                        m1 = size(U, 1);
                        n1 = size(U, 2);
    
                        U_11_12 = U(1:b1, :);
                        [~, ~, P_lupp, flops_lu] = fastLU(U_11_12');
                        U = U * P_lupp';
                        flops_perm = 2 * m1 * n1^2;
                        flops = flops + flops_lu + flops_perm;
                    case 'random'
                        % randomly permute the rotational matrix, as a comparison to QRCP and LUPP
                        P = eye(size(U, 2));
                        P = P(randperm(size(U, 2)), :);
                        U = U * P;
                        flops = flops + 2 * size(U, 1)^3;
                    otherwise
                        error("Invalid diagonalization method. Choose 'eig', 'qrcp', 'lupp', 'fastLU', or 'random'.");
                end
    
                % Update A using block transformation
                U_full = eye(n);
                rows = [blockIdx_i, blockIdx_j];
                cols = [blockIdx_i, blockIdx_j];
                U_full(rows, cols) = U;
                A = U_full' * A * U_full;
                A = (A + A') / 2; % Ensure symmetry
                flops = flops + 2 * 2 * n * (length(blockIdx_i) + length(blockIdx_j))^2;
    
                % Update Q
                Q = Q * U_full;
                flops = flops + 2 * n * (length(blockIdx_i) + length(blockIdx_j))^2;

            end
        end

        off_diag_norm = normOffDiag(A);
        fro_norm = FroNormOffDiag(A);
        sweep_OffNorm_history = [sweep_OffNorm_history; flops, sweeps, off_diag_norm, fro_norm];
        
        if off_diag_norm < eps_threshold
            break;
        end
    end

    % Return Diagonal Matrix and Orthonormal Matrix
    D = A;
end

function blockPairs = generateRowCyclicPairs(m)
    % Generate row-cyclic ordering
    blockPairs = [];
    for i = 1:m-1
        for j = i+1:m
            blockPairs = [blockPairs; i, j];
        end
    end
end

function blockPairs = generateColumnCyclicPairs(m)
    % Generate column-cyclic ordering
    blockPairs = [];
    for j = 2:m
        for i = 1:j-1
            blockPairs = [blockPairs; i, j];
        end
    end
end

function blockPairs = generateRandomPairs(m)
    % Generate all valid block pairs (i, j) with 1 <= i < j <= m
    pairs = [];
    for i = 1:m-1
        for j = i+1:m
            pairs = [pairs; i, j];
        end
    end

    % Shuffle the order of pairs randomly
    numPairs = size(pairs, 1);
    randomOrder = randperm(numPairs);
    blockPairs = pairs(randomOrder, :);
end
