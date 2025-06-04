function [U,Sigma,V] = blockOneSidedJacobi(G, b, tol, maxSweeps, wantVectors)
% BLOCKONESIDEDJACOBI_SVD  Approximate SVD of G via one-sided block Jacobi (Algorithm 5).
%
%   [U,Sigma,V] = blockOneSidedJacobi_SVD(G,b,tol,maxSweeps,wantVectors)
%
% INPUTS:
%   G          - n x n real matrix to be factorized (assumed square for this example).
%   b          - block size, must divide n exactly.
%   tol        - stopping tolerance for "near orthogonality".
%   maxSweeps  - maximum number of outer sweeps.
%   wantVectors- if true, accumulates rotations in V and outputs U,Sigma,V. 
%                otherwise, just sets V=[] and returns U=[],Sigma from G's column norms.
%
% OUTPUTS:
%   U     - approx left singular vectors. (Computed as G * inv(Sigma) at the end)
%   Sigma - diagonal matrix of singular values (column norms of G).
%   V     - approx right singular vectors.  (Only if wantVectors==true)
%
% EXAMPLE:
%   A = randn(8);  [U,S,V] = blockOneSidedJacobi_SVD(A,2,1e-12,20,true);

%------------------------------------------------------
% (1) INITIALIZATION
%------------------------------------------------------
[n,n2] = size(G);                                 % G is n x n
assert(n==n2,'Must be square for this example.');
assert(mod(n,b)==0,'Block size b must divide n exactly.');
nb = n / b;                                       % Number of blocks

% Initialize V if we want singular vectors
if wantVectors
    V = eye(n);
else
    V = [];
end

%------------------------------------------------------
% (2) OUTER LOOP of SWEEPS
%------------------------------------------------------
for sweep = 1:maxSweeps
    
    % (2a) Check if columns of G are "nearly orthonormal":
    Gram = G'*G;                                  % current Gram matrix
    offDiagFrob = norm(triu(Gram,1),'fro');       % measure off-diagonal
    if offDiagFrob < tol, break; end              % if small enough, done
    
    changed = false;                              % track if any rotation applied
    
    % (2b) FOR all block pairs 1 <= I < J <= nb
    for I = 1:(nb-1)
        for J = I+1:nb
            
            % Indices for block columns [ (I-1)*b+1 : I*b ] etc.
            colsI = (I-1)*b + (1:b);
            colsJ = (J-1)*b + (1:b);
            
            % (3) FORM the 2b x 2b submatrix A_hat = G(:, [I,J])^T G(:, [I,J]).
            %     We do this "on the fly":
            Gsub = [G(:, colsI), G(:, colsJ)];    % shape is n x 2b
            A_hat = Gsub' * Gsub;                % shape is 2b x 2b
            
            % (4) Decide if A_hat is "far from diagonal": 
            offSub = norm(triu(A_hat,1),'fro');
            diagSub= norm(diag(A_hat),2);
            if offSub > tol*diagSub   % "far enough" => do a rotation
                
                % (5) Compute eig of A_hat, i.e. A_hat ~ V_hat * D_hat * V_hat'
                %     (for a more accurate SVD, do 'svd' on A_hat)
                A_hat = (A_hat + A_hat')/2;      % ensure symmetric
                [V_hat, ~] = eig(A_hat);        % ~ for eigenvalues
                
                % (optional) ensure det(V_hat) > 0 for a rotation
                if det(V_hat) < 0
                    V_hat(:,1) = -V_hat(:,1);
                end
                
                % (6) Apply that 2b x 2b rotation to block columns I,J of G:
                newGsub = Gsub * V_hat;         % shape n x 2b
                G(:, colsI) = newGsub(:,1:b);
                G(:, colsJ) = newGsub(:,b+1:end);
                
                % (7) If singular vectors requested, also update V:
                if wantVectors
                    Vsub = [V(:, colsI), V(:, colsJ)] * V_hat; 
                    V(:, colsI) = Vsub(:,1:b);
                    V(:, colsJ) = Vsub(:,b+1:end);
                end
                
                changed = true;
            end
        end
    end
    
    % If no block changed, we can stop early
    if ~changed, break; end
end

%------------------------------------------------------
% (8) BUILD final Sigma from column norms, and U from G
%------------------------------------------------------
% Sigma is the diagonal matrix of column norms of G
colNorms = zeros(n,1);
for j=1:n
    colNorms(j) = norm(G(:,j),2);
end
Sigma = diag(colNorms);

% If we want left singular vectors: U = G * inv(Sigma)
if wantVectors
    U = G;
    for j=1:n
        if colNorms(j) > 0
            U(:,j) = U(:,j) / colNorms(j);
        end
    end
else
    U = [];
end
end
