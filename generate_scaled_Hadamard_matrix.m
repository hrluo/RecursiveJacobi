function A = generate_scaled_Hadamard_matrix(n)
% Generate matrix A with A = Q*D*Q', 
% where D is a random diagonal matrix and Q is a scaled Hadamard matrix (i.e., Q = (1/sqrt(n))*hadamard(n)).
%  
% Input: matrix size n
% Output: matrix A

% Check input validity
isPowerOfTwo = @(x) x > 0 && mod(x, 1) == 0 && (bitand(x, x-1) == 0);

% check whether n, n/12, or n/20 is a power of 2
if isPowerOfTwo(n) || isPowerOfTwo(n/12) || isPowerOfTwo(n/20)
    disp('n satisfies the condition.');
else
    error('The variable n must be an integer, and n, n/12, or n/20 must be a power of 2.');
end

diag_elements = randn(1, n);
D = diag(diag_elements);
Q = (1/sqrt(n))*hadamard(n);

A = Q * D * Q';
% To ensure output A is symmetric
A = (A + A')/2;
end