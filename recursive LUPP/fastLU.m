function [L, U, P, flops] = fastLU(A, flops)
    if nargin < 2
        flops = 0;  % Initialize flops if not provided
    end

    [n, m] = size(A);  % Get the dimensions of matrix A
    
    % Ensure the matrix is square or tall-and-skinny
    if n < m
        error('Matrix must be square or tall-and-skinny');
    end

    % Initialize permutation matrix as identity
    P = eye(n);

    % Base case: if matrix is a single column
    if m == 1
        [~,pivot] = max(abs(A));
        if pivot ~= 1
            % Swap first element with pivot (count swap flops)
            temp = A(1);
            A(1) = A(pivot);
            A(pivot) = temp;
            % For permutation matrix update, we might count some basic operations
            % But for simplicity, we can ignore small overhead here
            P(1,1) = 0;  P(1,pivot) = 1;
            P(pivot,1) = 1;  P(pivot,pivot) = 0;
            flops = flops + n;  % Rough cost for finding and swapping
        end
        % Base calculations for L and U
        L = A / A(1);
        U = A(1);
        % Small cost for division, assignment, etc...
        flops = flops + 5;  % arbitrary small constant for operations
    else
        % Partial Pivoting: Find the pivot row in the first column
        [~, pivot] = max(abs(A(:,1)));
        if pivot ~= 1
            % Swap rows in A
            A([1 pivot], :) = A([pivot 1], :);
            % Update permutation matrix P
            P([1 pivot], :) = P([pivot 1], :);
            flops = flops + 2*m;  % cost for swapping rows of size m
        end

        % Recursive case: split matrix into submatrices
        mid = floor(m / 2);
        midU = ceil(m / 2);

        % Step (a): Perform LU on the left half block
        [L_L, U_L, P_L, flops] = fastLU(A(1:n, 1:mid), flops);

        % Update overall permutation matrix
        P_temp = eye(n);
        P_temp(1:n,1:n) = P_L;
        P = P_temp * P;
        flops = flops + n^2;  % cost for applying permutation matrix

        % Multiply back half of A by P_L'
        A(1:n,(mid+1):m) = P_L'*A(1:n,(mid+1):m);
        flops = flops + 2*n*(m-mid);  % cost for multiplication by 2x2-like but generalized as dense multiply

        % Step (b): Update the upper right corner of A
        % Invert the leading L_L block and multiply
        L_inv = inv(L_L(1:mid,1:mid));
        flops = flops + (2/3)*(mid^3);  % cost for inversion
        A(1:mid, (mid+1):m) = L_inv * A(1:mid, (mid+1):m);
        flops = flops + 2*(mid^2)*(m-mid);  % multiply cost

        % Step (c): Update the Schur Complement
        A((mid+1):n,(mid+1):m) = A((mid+1):n,(mid+1):m) ...
            - L_L((mid+1):n,1:mid)*A(1:mid,(mid+1):m);
        flops = flops + 2*(n-mid)*mid*(m-mid);  % cost for this rank-update

        % Step (d): Perform LU on Schur complement
        [L_R, U_R, P_R, flops] = fastLU(A((mid+1):n, (mid+1):m), flops);

        % Update permutation matrix for Schur complement
        P_temp = eye(n);
        P_temp((mid+1):n, (mid+1):n) = P_R;
        P = P_temp * P;
        flops = flops + (n-mid)^2;  % cost for applying permutation

        % Step (e-f): Assemble L and U matrices
        L = horzcat(P_temp*L_L, [zeros(mid, midU); L_R]);
        U = vertcat(horzcat(U_L, A(1:mid,(mid+1):m)),horzcat(zeros(midU,mid), U_R));
        % Costs for assembly are relatively minor compared to major steps
        flops = flops + n*m;  % rough overhead for concatenation/assignment
    end
end

