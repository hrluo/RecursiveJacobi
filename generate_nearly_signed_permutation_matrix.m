function A = generate_nearly_signed_permutation_matrix(n, delta)
% Generate matrix A with A = Q*D*Q', 
% where D is a random diagonal matrix and Q is nearly a signed permutation matrix (i.e., [Q,~] = qr(eye(n) + delta*randn(n)) for some small delta > 0)
%  
% Input: matrix size n, perturbation parameter delta
% Output: matrix A

diag_elements = randn(1, n);
D = diag(diag_elements);
[Q, ~] = qr(eye(n) + delta * randn(n));

A = Q * D * Q';
% To ensure output A is symmetric
A = (A + A')/2;
end