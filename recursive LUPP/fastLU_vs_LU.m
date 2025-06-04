% Define the range of matrix sizes
ns = 128:128:2048;  % adjust as needed for performance measurement
fastLU_flops = zeros(size(ns));
standardLU_flops = zeros(size(ns)); 

% Loop over each matrix size
for k = 1:length(ns)
    n = ns(k);
    A = rand(n);
    
    % Use custom fastLU function that returns flop count
    [L, U, P, totalFlops] = fastLU(A);
    fastLU_flops(k) = totalFlops;
    
    % Theoretical flop count for standard LU without pivoting
    standardLU_flops(k) = (2/3) * n^3;%(1/6)*n*(4*n+1)*(n-1)
end

% Plotting the results
figure;
plot(ns, fastLU_flops, 'b-o', 'LineWidth', 2); hold on;
plot(ns, standardLU_flops, 'r--s', 'LineWidth', 2);
xlabel('Matrix Size n');
ylabel('Flop Count');
legend('fastLU Flop Count', 'Theoretical Standard LU Flop Count (2/3 n^3)', ...
       'Location', 'NorthWest');
title('Comparison of Flop Counts for LU Decompositions');
grid on;
