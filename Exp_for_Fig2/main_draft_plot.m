clear

seed = 101;
rng(seed);

% Example usage
n_size = 512;
num_blocks_values = [4, 8, 16, 32];  % Different num_blocks values for Block Jacobi
f_values = [0.2, 0.4, 0.6, 0.8];  % Different f values for Recursive Jacobi

% Randomly generated symmetric matrix
A = randn(n_size, n_size);
A = (A + A') / 2;  % Make A symmetric

% Set parameters
n_threshold = 4;
eps_threshold = 1e-7 * max(max(abs(A)));
n3_ratio = Inf;   % Ensure convergence for every method

% Record data
results_classical = struct();
results_block = struct();
results_recursive = struct();

%% Classical Jacobi
fprintf("Classical Jacobi\n");
tic;
A1 = A;
[Q_classical, D_classical, flops_classical, sweeps_classical, sweep_history_classical] = ...
    classicalJacobi(A1, eps_threshold, 'trig', n3_ratio);
time_classical = toc;

% record the data
results_classical.Q = Q_classical;
results_classical.D = D_classical;
results_classical.flops = flops_classical;
results_classical.sweeps = sweeps_classical;
results_classical.sweeps_history = sweep_history_classical;
results_classical.time = time_classical;

% block Jacobi
for i = 1:length(num_blocks_values)
    num_blocks = num_blocks_values(i);
    blockSizes = repmat(n_size / num_blocks, 1, num_blocks);
    fprintf("Block Jacobi with num_blocks=%d\n", num_blocks);
    
    tic;
    A1 = A;
    [Q_block, D_block, flops_block, sweeps_block, sweep_history_block] = ...
        BlockJacobi(A1, blockSizes, eps_threshold, "rowcyclic", "eig", n3_ratio);
    time_block = toc;
    
    % record the data
    results_block(i).numblocks = num_blocks_values(i);
    results_block(i).Q_block = Q_block;
    results_block(i).D_block = D_block;
    results_block(i).flops_block = flops_block;
    results_block(i).sweeps_block = sweeps_block;
    results_block(i).sweep_history_block = sweep_history_block;
    results_block(i).time = time_block;

end


% Recursive Jacobi - Varying f
for i = 1:length(f_values)
    f_val = f_values(i);
    fprintf("Recursive Jacobi with f=%0.1f\n", f_val);
    
    tic;
    A1 = A;
    [Q_recursive, D_recursive, flops_recursive, sweeps_recursive, sweep_history_recursive, ~] = ...
        RecursiveJacobiplain(A1, n_threshold, f_val, eps_threshold, 0, 0, n3_ratio); 
    time_recursive = toc;
    
    %  record the data
    results_recursive(i).f = f_val;
    results_recursive(i).Q_recursive = Q_recursive;
    results_recursive(i).D_recursive = D_recursive;
    results_recursive(i).flops_recursive = flops_recursive;
    results_recursive(i).sweeps_recursive = sweeps_recursive;
    results_recursive(i).sweep_history_recursive = sweep_history_recursive;
    results_recursive(i).time = time_recursive;

end

%% Plotting
% === Set layout dimensions ===
W_tile = 1200;   % Width of each tile
H_tile = 600;   % Height of each tile
layout = [2, 1];
num_rows = layout(1);
num_cols = layout(2);
fig_width = W_tile * num_cols + 10;
fig_height = H_tile * num_rows + 10;

% === Create figure ===
figure;
set(gcf, 'Position', [100, 0, fig_width, fig_height]);
tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

% === Set colors and linestyles ===
line_styles = {'-', '--', ':', '-.'};
colors_block = lines(length(results_block));
colors_recursive = lines(length(results_recursive));
markers = {'*', 's', 'd', '^'};

% === Norm types ===
norm_types = {'maxabs', 'fro'};  % column 3 and 4 in sweep_history

% ---------------------- Row 1: Classical + Block ----------------------
for col = 1:1
    nexttile(col);

    % Plot Classical Jacobi
    sweeps_history = results_classical.sweeps_history;
    unique_sweeps = sweeps_history(:, 2);
    
    flop_data = sweeps_history(:, 1);       % extract flops
    error_data = sweeps_history(:, 3);      % extract max abs off-diagonal value

    loglog(flop_data, error_data, ...
                'LineStyle', '--', ...
                'Marker', 'o', ...
                'Color', [0 0 0], ...
                'MarkerFaceColor', [0 0 0], ...
                'LineWidth', 2, ...
                'DisplayName', 'Scalar Jacobi');

    hold on;

    % Annotate Classical
    text(0.97 * flop_data(end), 0.3 * error_data(end), num2str(unique_sweeps(end)), ...
        'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');

    % Plot Block Jacobi
    for i = 1:length(results_block)
        
        % Plot Classical Jacobi
        sweeps_history = results_block(i).sweep_history_block;
        unique_sweeps = sweeps_history(:, 2);
        
        flop_data = sweeps_history(:, 1);       % extract flops
        error_data = sweeps_history(:, 3);      % extract max abs off-diagonal value

        loglog(flop_data, error_data, ...
            'LineStyle', ':', ...
            'Marker', markers{i}, ...
            'Color', colors_block(i,:), ...
            'MarkerFaceColor', colors_block(i,:), ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('Block Jacobi (%d blocks)', results_block(i).numblocks));

        if i == 1
            text(flop_data(end), 0.45 * error_data(end), ...
                 num2str(unique_sweeps(end)), 'Color', colors_block(i,:), ...
                 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        elseif i == 2
            text(1.05 * flop_data(end), 0.4 * error_data(end), ...
                 num2str(unique_sweeps(end)), 'Color', colors_block(i,:), ...
                 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        elseif i == 3
            text(1.05 * flop_data(end), 0.4 * error_data(end), ...
                 num2str(unique_sweeps(end)), 'Color', colors_block(i,:), ...
                 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        elseif i == 4
        text(1.1 * flop_data(end), 0.4 * error_data(end), ...
             num2str(unique_sweeps(end)), 'Color', colors_block(i,:), ...
             'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end

    % Decorations
    xlabel('FLOPs');

    if col == 1
        ylabel('Max Abs Off-Diagonal');
    else
        ylabel('Frobenius Off-Diagonal');
    end

    title(sprintf('Scalar & Block Jacobi (%s)', norm_types{col}));
    
    if col == 1
        yline(eps_threshold, '--r', ...
            'LineWidth', 1.5, ...
            'DisplayName', '1e-7 × max(max(abs(A)))');
    else
        eps_norm_threshold = 1e-7 * norm(A, "fro");
        yline(eps_norm_threshold, '--r', '1e-7 * norm(A, "fro")', ...
            'LabelHorizontalAlignment','right', 'LineWidth', 1.5, ...
            'HandleVisibility','off');
    end

    grid on;
    legend('Location', 'northeast'); legend('boxoff'); legend('AutoUpdate','off');
    set(legend, ...
    'EdgeColor', 'black', ...   % edge color
    'Color', 'white', ...       % background color
    'Box', 'on');               % display the edge

    xlim([7e8, 4e11]);
    ylim([1e-8, 1e1]);
end

% ---------------------- Row 2: Recursive Jacobi ----------------------
% for col = 1:2
%     nexttile(2 + col);  % Tiles 3 and 4
for col = 1:1
    nexttile(1 + col);  % Tiles 3 and 4

    for i = 1:length(results_recursive)
        
        % Plot Classical Jacobi
        sweeps_history = results_recursive(i).sweep_history_recursive;
        unique_sweeps = sweeps_history(:, 2);
        
        flop_data = sweeps_history(:, 1);       % extract flops
        error_data = sweeps_history(:, 3);      % extract max abs off-diagonal value

        loglog(flop_data, error_data, ...
                    'LineStyle', '-', ...
                    'Marker', markers{i}, ...
                    'MarkerFaceColor', colors_recursive(i,:), ...
                    'Color', colors_recursive(i,:), ...
                    'LineWidth', 2, ...
                    'DisplayName', sprintf('f = %.1f', results_recursive(i).f));
        
        hold on;
        
        text(flop_data(end), 0.4 * error_data(end), ...
             num2str(unique_sweeps(end)), 'Color', colors_recursive(i,:), ...
             'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

    % Decorations
    xlabel('FLOPs');
    
    if col == 1
        ylabel('Max Abs Off-Diagonal');
    else
        ylabel('Frobenius Off-Diagonal');
    end

    title(sprintf('Recursive Jacobi (%s)', norm_types{col}));
    
    if col == 1
        yline(eps_threshold, '--r', ...
            'LineWidth', 1.5, ...
            'DisplayName', '1e-7 × max(max(abs(A)))');
    else
        eps_norm_threshold = 1e-7 * norm(A, "fro");
        yline(eps_norm_threshold, '--r', '1e-7 * norm(A, "fro")', ...
            'LabelHorizontalAlignment','right', 'LineWidth', 1.5, ...
            'HandleVisibility','off');
    end

    grid on;
    legend('Location', 'southwest'); legend('boxoff'); legend('AutoUpdate','off');
    set(legend, ...
    'EdgeColor', 'black', ...   % edge color
    'Color', 'white', ...       % background color
    'Box', 'on');               % display the edge
    
    xlim([7e8, 4e11]);
    ylim([1e-8, 1e1]);
end
