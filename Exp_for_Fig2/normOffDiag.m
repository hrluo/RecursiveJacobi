function norm_off_diag = normOffDiag(A)
% Computes the maximum off-diagonal absolute value of a symmetric matrix A
    norm_off_diag = max(max(abs(triu(A, 1))));
end
