function [Q, D, flops, sweeps, sweep_OffNorm_history, new_break_flag] = ...
    RecursiveJacobiplain(A, n_threshold, f, eps_threshold, recdepth, break_flag, ...
                              n3_ratio)
    % vanilla Recursive Jacobi Algorithm for the Symmetric Eigenproblem
    % Inputs:
    %   A - symmetric matrix
    %   n_threshold - size threshold to switch to direct eigenvalue calculation
    %   f - log block size
    %   eps_threshold - tolerance for stopping condition
    %   recdepth - current recursion depth
    %   break_flag - flag for termination
    %   n3_ratio - threshold ratio for flops
    % Outputs:
    %   Q - orthogonal matrix of eigenvectors
    %   D - near-diagonal matrix of eigenvalues
    %   flops - total FLOPs performed
    %   sweeps - total sweeps performed
    %   sweep_OffNorm_history - a log of execution, recorded in the format [flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm], with entries captured at the end of each sweep
    %   new_break_flag - flag for termination


    % Ensure A is a symmetric matrix
    if break_flag == 1
        new_break_flag = 1;
        return;
    end
    new_break_flag = 0;
    A = (A + A') / 2; % Make sure A is symmetric
    D = A;
    assert(issymmetric(A), 'Matrix A must be symmetric');

    % Get the size of the matrix A
    n = size(A, 1);

    % Initialize Q as the identity matrix
    Q = eye(n);

    % Initialize FLOPs
    flops = 0;

    % Initialize sweeps
    sweeps = 0;

    % Initialize recording histories
    sweep_OffNorm_history = [];

    % Record if recdepth == 0
    if recdepth == 0
        sweep_OffNorm_history = [sweep_OffNorm_history; 0, 0, normOffDiag(A), FroNormOffDiag(A)];
    end
    
    % Base case: if n is small enough, solve eigenproblem directly
    if n <= n_threshold
        [Q, D] = eig(A);
        D = diag(diag(D));
        flops = flops + (8 + 2/3) * n^3;

        % If recdepth == 0, update the log history
        if recdepth == 0
    	    sweeps = sweeps + 1;
            sweep_OffNorm_history = [sweep_OffNorm_history; flops, sweeps, normOffDiag(Q'*A*Q), FroNormOffDiag(Q'*A*Q)];           
        end
        
        return;
    end

    % Determine the blocking parameter b
    b = round(n^f);
    if b >= round(n / 2)
        [Q, D] = eig(A);
        D = diag(diag(D));
        flops = flops + (8 + 2/3) * n^3;

        % If recdepth == 0, update the log history
        if recdepth == 0
            sweeps = sweeps + 1;
            sweep_OffNorm_history = [sweep_OffNorm_history; flops, sweeps, normOffDiag(Q'*A*Q), FroNormOffDiag(Q'*A*Q)];
        end

        return;
    end 

    % Repeat until off-diagonal entries are small enough
    while flops < n3_ratio * n^3
        % Only count sweeps at top level
        if recdepth == 0
            sweeps = sweeps + 1;
        end
        for I = 1:b:(n-b+1)
            for J = (I+b):b:n

                % Define indices for current block within the submatrix
                I_block_indices = I:min(I+b-1, n);
                J_block_indices = J:min(J+b-1, n);
                rows = [I_block_indices, J_block_indices];
                cols = rows;  % symmetric matrix


                % Extract 2b-by-2b submatrix A_hat
                A_hat = A(rows, cols);
                A_hat = (A_hat + A_hat') / 2;

                if normOffDiag(A_hat) >= eps_threshold
                    % Recursive call with updated global indices
                    [Q_hat, ~, sub_flops, sub_sweeps, ~, break_flag1] = ...
                        RecursiveJacobiplain(A_hat, n_threshold, f, eps_threshold, ...
                                                  recdepth + 1, break_flag, ...
                                                  n3_ratio);
                    new_break_flag = break_flag1;
                    flops = flops + sub_flops;
    
                    % Update A and Q with block transformations
                    A(rows, :) = Q_hat' * A(rows, :);
                    A(:, cols) = A(:, cols) * Q_hat;
    
                    % Update the corresponding columns of Q
                    Q(:, rows) = Q(:, rows) * Q_hat;
    
                    % Update total flops
                    flops = flops + 3 * (2 * b)^2 * n;
                end
                
            end
        end

        % After finishing one sweep (at the top level), update log history once.
        if recdepth == 0
            sweep_OffNorm_history = [sweep_OffNorm_history; ...
                                     flops, sweeps, normOffDiag(A), FroNormOffDiag(A)];
        end
        
        % Update D correspondingly
        D = A;

        % Update convergence condition and off-diagonal elements
        max_off_diag = normOffDiag(A); 
        if max_off_diag < eps_threshold
            new_break_flag = 1;
            break;
        end
    end
end

