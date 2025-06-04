function [Q, D, flops, sweeps, sweep_OffNorm_history] = classicalJacobi(A, eps_threshold, method, n3_ratio)
    % Classical (scalar) Jacobi with row-cyclic ordering for the Symmetric Eigenproblem with multiple methods
    % Inputs:
    %   A - symmetric matrix
    %   eps_threshold - tolerance for stopping condition
    %   method - 'eig'              % use Matlab eig() to diagonalize the 2-by-2 sub-matrix
    %            'trig'             % diagonalize the 2-by-2 sub-matrix by explicitly calculating the rotation angle
    %            'adversarial'      % deliberately adding pi/2 to the rotation angle to resist convergence, 
    %                               % based on THE Cyclic Jacobi Method for Computing the Principal Values of a Complex Matrix, by Forsythe and Henrici
    %   n3_ratio - threshold for limiting the maximum FLOPs
    % Outputs:
    %   Q - orthogonal matrix of eigenvectors
    %   D - near-diagonal matrix of eigenvalues
    %   flops - total FLOPs performed
    %   sweeps - total sweeps performed
    %   sweep_OffNorm_history - a log of execution, recorded in the format [flops, sweeps, maximum off-diagonal absolute value, off-diagonal Frobenius norm], with entries captured at the end of each sweep

    % Ensure A is symmetric
    A = (A + A') / 2;
    assert(issymmetric(A), 'Matrix A must be symmetric');

    % Initialize parameters
    n = size(A, 1);
    Q = eye(n);
    flops = 0;
    sweeps = 0;
    sweep_OffNorm_history = [flops, sweeps, normOffDiag(A), FroNormOffDiag(A)];

    % Main iteration loop
    while flops < n3_ratio * n^3
        sweeps = sweeps + 1;
        
        for i = 1:n-1
            for j = i+1:n
                if abs(A(i, j)) >= eps_threshold
                    if strcmp(method, 'eig')
                        % Eigen decomposition
                        A_hat = A([i, j], [i, j]);
                        A_hat = (A_hat + A_hat') / 2;
                        [V, D_hat] = eig(A_hat);
                        Q_hat = V;

                        % Ensure Q_hat is a proper rotation matrix
                        tolerance = 1e-10;
                        if abs(det(Q_hat) - 1) < tolerance
                            R = Q_hat;
                        elseif abs(det(Q_hat) + 1) < tolerance
                            R = Q_hat;
                            R(:,1) = -R(:,1);  % Adjust to make determinant +1
                        else
                            error('Determinant of Q_hat is not ±1. Cannot proceed.');
                        end

                    elseif strcmp(method, 'trig')
                        % Trigonometric transformation
                        theta_rad = atan((2 * A(i, j)) / (A(i, i) - A(j, j))) / 2;
                        cos_theta = cos(theta_rad);
                        sin_theta = sin(theta_rad);
                        R = [cos_theta, sin_theta; -sin_theta, cos_theta];
                    elseif strcmp(method, 'adversarial')
                        % We still do trigonometric transformation, but
                        % we deliberately set the rotation angle to be
                        % out of any open intervals within [-pi/2, pi/2]
                        theta_rad = (atan((2 * A(i, j)) / (A(i, i) - A(j, j))) + pi) / 2;
                        cos_theta = cos(theta_rad);
                        sin_theta = sin(theta_rad);
                        R = [cos_theta, sin_theta; -sin_theta, cos_theta];
                    else
                        error("Unknown method. Use 'eig', 'trig', or 'adversarial'.");
                    end

                    % Apply rotation to A and Q for indices (i, j)
                    A([i, j], :) = R * A([i, j], :);
                    A(:, [i, j]) = A(:, [i, j]) * R';
                    Q(:, [i, j]) = Q(:, [i, j]) * R';

                    % Symmetrize A to ensure convergence
                    A = (A + A') / 2;

                    % Update FLOPs
                    flops = flops + 2 * 2^3 + 3 * 2 * 3 * n;    % flops = flops + flops of eigen-decomposition + flops of updating the whole matrix
                end
            end
        end

        % Update sweep_OffNorm_history
        off_diag_norm = normOffDiag(A);
        fro_norm = FroNormOffDiag(A);
        sweep_OffNorm_history = [sweep_OffNorm_history; flops, sweeps, off_diag_norm, fro_norm];

        if off_diag_norm < eps_threshold
            break;
        end
    end

    % Output the diagonalized matrix D
    D = A;
end
