function fro_norm_off_diag = FroNormOffDiag(A)
% Computes the Frobenius norm of the off-diagonal elements of a symmetric matrix A
    fro_norm_off_diag = sqrt(2*sum(sum(triu(A, 1).^2)));
end
