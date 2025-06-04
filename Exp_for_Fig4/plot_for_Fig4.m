%% Parameters and Initialization
clear

seed = 101;
rng(seed);

n_size = 512;
n_threshold = 4;
f = 0.4;
n3_ratio = 1000;

A = randn(n_size, n_size);
A = (A + A') / 2;

eps_threshold = 1e-7 * max(max(abs(A)));   % Convergence criteria
eps_norm_threshold = 1e-7 * norm(A, "fro");   % for Comparison
flops_limit = n3_ratio * n_size^3;

% diffeent bottom case precisions
bottom_case_precisions = [1e-5, 1e-7, 1e-9, 1e-11, 1e-13];

results_recursive = struct();

%% Run experiments for different bottom case precisions
for i = 1:length(bottom_case_precisions)
    bottom_case_precision = bottom_case_precisions(i);
    bottom_case_precision = bottom_case_precision * max(max(abs(A)));

    fprintf("Running with bottom case precision = %e\n", bottom_case_precisions(i));
    
    A1 = A;
    tic;
    [Q, D, flops, sweeps, sweep_OffNorm_history, new_break_flag] = ...
        RecursiveJacobiplain(A1, n_threshold, f, eps_threshold, 0, 0, ...
                             n3_ratio, bottom_case_precision);
    time = toc;

    results_recursive(i).bottom_case_precision = bottom_case_precisions(i);
    results_recursive(i).Q = Q;
    results_recursive(i).D = D;
    results_recursive(i).flops = flops;
    results_recursive(i).sweeps = sweeps;
    results_recursive(i).history = sweep_OffNorm_history;
    results_recursive(i).time = time;
end

%% Plotting
% Fixed size for each subplot tile (in pixels)
W_tile = 600;   % Width of each tile
H_tile = 300;   % Height of each tile

layout = [1, 2];
% layout = [2, 1];
num_rows = layout(1);
num_cols = layout(2);

% Calculate total figure size based on tile size and layout
fig_width = W_tile * num_cols + 10;   % Add some margin
fig_height = H_tile * num_rows + 10;  % Add some margin

% Create figure window and set its size
figure;
set(gcf, 'Position', [100, 0, fig_width, fig_height]);
t = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

% Define norm types: column index offset (3 for maxabs, 4 for fro)
norm_types = {'maxabs', 'fro'};
line_styles = {'-', '--', ':', '-.'};  % Different line styles for each precision
colors = lines(4);  % Use distinct colors
markers = {'*', 's', 'd', '^'};

legend_handles = [];
legend_labels = {};

for norm_idx = 1:2
    nexttile;
    first = 0;

    for j = 1:4  % Loop over precision values
        % Extract sweep history
        sweep_history = results_recursive(j).history;
        unique_sweeps = sweep_history(:, 2);

        % Initialize data containers
        flop_data = sweep_history(:, 1);
        error_data = sweep_history(:, norm_idx + 2);   % Choose max abs and fro norm accordingly

        % Plot FLOPs vs. error in log-log scale
        % Format legend label
        precision_str = sprintf('%s', bottom_case_precisions(j));
        precision_str = strrep(precision_str, '0', '');
        precision_str = strrep(precision_str, '.', '');
        
        h = loglog(flop_data, error_data, 'Marker', markers{j}, ...
                'MarkerFaceColor', colors(j,:), ...
                'LineStyle', line_styles{j}, ...
                'Color', colors(j,:), ...
                'LineWidth', 2);
        
        if norm_idx == 1
            legend_handles(end+1) = h;
            legend_labels{end+1} = sprintf('precision=%s', precision_str);
        end

        if first == 0
            hold on;
            first = 1;
        end

        % Add label to the last sweep
        if j == 2
            text(0.9*flop_data(end), 0.3 * error_data(end), num2str(unique_sweeps(end)), ...
                'Color', colors(j,:), 'FontSize', 10, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        elseif j == 4
            text(1.1*flop_data(end), 0.3 * error_data(end), num2str(unique_sweeps(end)), ...
                'Color', colors(j,:), 'FontSize', 10, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        else
            text(flop_data(end), 0.3 * error_data(end), num2str(unique_sweeps(end)), ...
                'Color', colors(j,:), 'FontSize', 10, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end

    % Set labels and title for each subplot
    xlabel('FLOPs', 'FontSize', 12);
    if norm_idx == 1
        ylabel('Max Abs Off-Diagonal', 'FontSize', 12);
        title('Recursive Jacobi with different bottom case precisions (maxabs)');
    else
        ylabel('Frobenius Off-Diagonal', 'FontSize', 12);
        title('Recursive Jacobi with different bottom case precisions (fro)');
    end

    % Plot decorations
    grid on;

    % Add horizontal threshold line
    if norm_idx == 1
        h = yline(eps_threshold, '--r', ...
                'LineWidth', 1.5);
        legend_handles(end+1) = h;
        legend_labels{end+1} = '1e-7 × max(max(abs(A)))';
    else
        h = yline(eps_norm_threshold, '-r', ...
                'LineWidth', 1.5);
        legend_handles(end+1) = h;
        legend_labels{end+1} = "1e-7 × norm(A, ""fro"")";
    end

    % Add vertical flops limit line
    h = xline(flops_limit, '--', ...
                'Color', [0.5, 0.5, 0.5], ...
                'LineWidth', 1.5);
    
    if norm_idx == 2
        legend_handles(end+1) = h;
        legend_labels{end+1} = "flops limit";
    end

    % Set axis limits
    xlim([1e9, 2e11]);
    ylim([5*1e-10, 1e3]);

    if norm_idx == 2
        lgd = legend(legend_handles, legend_labels);
        lgd.Layout.Tile = 'east';
        lgd.Box = 'on';
    end

end

