n = 1024;
f = 0.8;
n_threshold = 4;
structure_array = Reveal_recursive_structure(n, n_threshold, f);
disp(structure_array);
depths = size(structure_array, 1); % Recursion depth

%% ---- 3D Visualization ----
figure;
hold on;

% Draw base level (level 0) n x n matrix
X = [0, n, n, 0];
Y = [0, 0, n, n];
Z = zeros(1,4);
fill3(X, Y, Z, 'w', 'EdgeColor', 'k', 'LineWidth', 2); % White fill, opaque

for i = 1:depths
    d = structure_array(i, 1); % Recursion level
    n_size = structure_array(i, 2); % Matrix size
    block_size = structure_array(i, 3); % Block size
    
    % Draw the outer boundary of the matrix (using fill3 to ensure opacity)
    X = [0, n_size, n_size, 0];
    Y = [0, 0, n_size, n_size];
    Z = d * ones(1,4);
    fill3(X, Y, Z, 'w', 'EdgeColor', 'k', 'LineWidth', 2); % White fill, opaque
    
    % Draw block boundaries only if recursion level is greater than 0
    if d > 0
        for j = block_size:block_size:n_size
            % Draw horizontal lines
            line([0, n_size], [j, j], [d, d], 'Color', 'k', 'LineWidth', 2);
            % Draw vertical lines
            line([j, j], [0, n_size], [d, d], 'Color', 'k', 'LineWidth', 2);
        end
    end
end

hold off;
view(45, 30); % Adjust 3D view
xlabel('matrix size');
ylabel('matrix size');
zlabel('Recursion Depth');
title('3D Recursive Block Structure');
grid on;

% Display only integer values on the Z-axis
zticks(1:depths);

