%% Short script to compare blockOneSidedJacobi to MATLAB's svd

clear; clc;

seed = 101;
rng(seed);

n = 64;                % matrix dimension
A = randn(n,n);       % random matrix
% Optionally make it ill-conditioned or rectangular, e.g. A = randn(20,6).

epsTol = 1e-16; 
blockSize = 16;
maxsweeps = 2000;  % Just a placeholder
n3_ratio = 100;

% Call one-sided Jacobi SVD function
A1 = A;
[U,Sigma,V] = blockOneSidedJacobi(A1, blockSize, epsTol, maxsweeps, true);

% Compare singular values to MATLAB's built-in SVD:
A1 = A;
[Utrue, Strue, Vtrue] = svd(A1);

% The diagonal entries of Strue are the "true" singular values. 
% The diagonal of Sigma is the approximate singular values. 
% Compare them in sorted order:
singValsApprox = diag(Sigma);        % from blockOneSidedJacobi
singValsTrue   = diag(Strue);

% Sort them if needed (some algorithms might produce them in any order).
[~, idxA] = sort(singValsApprox,'descend');
[~, idxT] = sort(singValsTrue,'descend');
svApproxSorted = singValsApprox(idxA);
svTrueSorted   = singValsTrue(idxT);
singValsApprox
singValsTrue
% Print or measure the difference in singular values
fprintf('Approx vs. True singular values:\n');
for k = 1:n
    fprintf('  sigma_approx = %9.5g   |   sigma_true = %9.5g   |   rel.error = %g\n',...
        svApproxSorted(k), svTrueSorted(k), ...
        abs(svApproxSorted(k) - svTrueSorted(k))/max(eps, svTrueSorted(k)));
end

