% Define the matrix sizes and num_blocks values
clear;

seed = 101;
rng(seed);

n_sizes = [128, 256, 512, 1024, 2048];  % Matrix sizes
num_blocks_values = [4, 8, 16, 32];  % num_blocks values
n3_ratio = Inf;  % Example ratio (adjust based on your requirements)


% Preallocate result struct array
% matrix_data(length(n_sizes)) = struct();
results(length(n_sizes), length(num_blocks_values)) = struct();

%% Running block Jacobi with different parameters
for i = 1:length(n_sizes)
    n_size = n_sizes(i);  % Current matrix size
    
    % Generate symmetric matrix
    A = randn(n_size, n_size);
    A = (A + A')/2;
    % matrix_data(i).matrix = A;

    % Set threshold
    eps_threshold = 1e-7 * max(max(abs(A)));

    for j = 1:length(num_blocks_values)
        num_blocks = num_blocks_values(j);  % Current number of blocks

        % Create equal-sized blocks
        blockSizes = repmat(n_size / num_blocks, 1, num_blocks);

        fprintf("Running Block Jacobi with num_blocks=%d and n_size=%d\n", num_blocks, n_size);

        A1 = A;
        % Run the Block Jacobi algorithm
        [Q_block, D_block, flops_block, sweeps_block, sweep_history_block] = ...
            BlockJacobi(A1, blockSizes, eps_threshold, "rowcyclic", "eig", n3_ratio);

        % Store everything into a struct (Here we only store sweeps for saving memory space)
        results(i, j).n_size = n_size;
        results(i, j).num_blocks = num_blocks;
        % results(i, j).A = A1;
        % results(i, j).eps_threshold = eps_threshold;
        % results(i, j).Q = Q_block;
        % results(i, j).D = D_block;
        % results(i, j).flops = flops_block;
        results(i, j).sweeps = sweeps_block;
        disp(sweeps_block); fprintf("\n");
        % results(i, j).sweep_history = sweep_history_block;
    end
end


%% Plotting
figure;
set(gcf, 'Position', [100, 100, 1600, 800]);  % Optional: wider figure

% Define line styles and colors
line_styles = {'-', '--', ':', '-.'};
colors = lines(length(num_blocks_values));

% Loop over each num_blocks value (fixed j)
for j = 1:length(num_blocks_values)
    sweep_vals = zeros(1, length(n_sizes));  % Preallocate sweep counts

    for i = 1:length(n_sizes)
        sweep_vals(i) = results(i, j).sweeps;
    end

    % Plot sweeps vs. n_sizes
    semilogx(n_sizes, sweep_vals, ...
        'LineStyle', line_styles{mod(j-1,length(line_styles))+1}, ...
        'Color', colors(j,:), ...
        'Marker', 'o', ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('num\\_blocks = %d', num_blocks_values(j)));
    
    hold on;

end

xlim([100,3000]);
ylim([1,10]);

% Get y-axis limits
y_limits = ylim;

for i = 1:length(n_sizes)
    % Plot vertical line without adding it to legend
    xl = xline(n_sizes(i), '--k', 'LineWidth', 0.8);
    xl.Annotation.LegendInformation.IconDisplayStyle = 'off';  % Hide from legend

    % Add label below the x-axis, slightly offset
    text(1.05*n_sizes(i), y_limits(1) + 0.3, sprintf('%d', n_sizes(i)), ...
        'FontSize', 11, ...                        % Larger font
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');              % Appear below axis
end


xlabel('Matrix Size (n)', 'FontSize', 12);
ylabel('Number of Sweeps', 'FontSize', 12);
title('Sweeps vs. Matrix Size for Block Jacobi', 'FontSize', 14);
legend('Location', 'northeast');
grid on;
hold off;
