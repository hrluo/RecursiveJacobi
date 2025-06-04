%% Parameters and initialization
clear

seed = 101;
rng(seed);

% Example usage
n_size = 512;
num_blocks = 256;  % Different num_blocks values for Block Jacobi
blockSizes = repmat(n_size / num_blocks, 1, num_blocks);

% Randomly generated symmetric matrix
A = randn(n_size, n_size);
A = (A + A') / 2;  % Make A symmetric

% Set parameters
eps_threshold = 1e-7 * max(max(abs(A)));
n3_ratio = 350;
flops_limit = n3_ratio * n_size^3;

% Run 5 Block Jacobi variants
results_blockjacobi = struct();

%%
% 1. Non-adversarial
A1 = A;
fprintf("Non-adversarial\n");
tic;
[Q_Nadv, D_Nadv, flops_Nadv, sweeps_Nadv, sweep_OffNorm_history_Nadv] = BlockJacobi_Nadv(A1, blockSizes, eps_threshold, "rowcyclic", "eig", n3_ratio);
time_Nadv = toc;
results_blockjacobi(1).label = 'Nadv';
results_blockjacobi(1).Q = Q_Nadv;
results_blockjacobi(1).D = D_Nadv;
results_blockjacobi(1).flops = flops_Nadv;
results_blockjacobi(1).sweeps = sweeps_Nadv;
results_blockjacobi(1).history = sweep_OffNorm_history_Nadv;
results_blockjacobi(1).time = time_Nadv;

% 2. Adversarial eig
A1 = A;
fprintf("Adversarial\n");
tic;
[Q_adv, D_adv, flops_adv, sweeps_adv, sweep_OffNorm_history_adv] = BlockJacobi_adv(A1, blockSizes, eps_threshold, "rowcyclic", "eig", n3_ratio);
time_adv = toc;
results_blockjacobi(2).label = 'Adv';
results_blockjacobi(2).Q = Q_adv;
results_blockjacobi(2).D = D_adv;
results_blockjacobi(2).flops = flops_adv;
results_blockjacobi(2).sweeps = sweeps_adv;
results_blockjacobi(2).history = sweep_OffNorm_history_adv;
results_blockjacobi(2).time = time_adv;

% 3. Adversarial QRCP
A1 = A;
fprintf("Adversarial QRCP\n");
tic;
[Q_adv_qrcp, D_adv_qrcp, flops_adv_qrcp, sweeps_adv_qrcp, sweep_OffNorm_history_adv_qrcp] = BlockJacobi_adv(A1, blockSizes, eps_threshold, "rowcyclic", "qrcp", n3_ratio);
time_adv_qrcp = toc;
results_blockjacobi(3).label = 'Adv-QRCP';
results_blockjacobi(3).Q = Q_adv_qrcp;
results_blockjacobi(3).D = D_adv_qrcp;
results_blockjacobi(3).flops = flops_adv_qrcp;
results_blockjacobi(3).sweeps = sweeps_adv_qrcp;
results_blockjacobi(3).history = sweep_OffNorm_history_adv_qrcp;
results_blockjacobi(3).time = time_adv_qrcp;

% 4. Adversarial LUPP
A1 = A;
fprintf("Adversarial LUPP\n");
tic;
[Q_adv_lupp, D_adv_lupp, flops_adv_lupp, sweeps_adv_lupp, sweep_OffNorm_history_adv_lupp] = BlockJacobi_adv(A1, blockSizes, eps_threshold, "rowcyclic", "lupp", n3_ratio);
time_adv_lupp = toc;
results_blockjacobi(4).label = 'Adv-LUPP';
results_blockjacobi(4).Q = Q_adv_lupp;
results_blockjacobi(4).D = D_adv_lupp;
results_blockjacobi(4).flops = flops_adv_lupp;
results_blockjacobi(4).sweeps = sweeps_adv_lupp;
results_blockjacobi(4).history = sweep_OffNorm_history_adv_lupp;
results_blockjacobi(4).time = time_adv_lupp;

% 5. Adversarial Random
A1 = A;
fprintf("Adversarial random\n");
tic;
[Q_adv_random, D_adv_random, flops_adv_random, sweeps_adv_random, sweep_OffNorm_history_adv_random] = BlockJacobi_adv(A1, blockSizes, eps_threshold, "random", "eig", n3_ratio);
time_adv_random = toc;
results_blockjacobi(5).label = 'Adv-Random';
results_blockjacobi(5).Q = Q_adv_random;
results_blockjacobi(5).D = D_adv_random;
results_blockjacobi(5).flops = flops_adv_random;
results_blockjacobi(5).sweeps = sweeps_adv_random;
results_blockjacobi(5).history = sweep_OffNorm_history_adv_random;
results_blockjacobi(5).time = time_adv_random;

%% Plotting
% plot size
fig_width = 1210;   % Width of the figure window in pixels
fig_height = 600;   % Height of the figure window in pixels
margin_px = 75;     % Desired margin around the plot (in pixels)

figure;
set(gcf, 'Position', [100, 100, fig_width, fig_height]);

% Convert pixel margin to normalized units (range from 0 to 1)
margin_x = margin_px / fig_width;
margin_y = margin_px / fig_height;

% Set the axes position: [left, bottom, width, height] in normalized units
% This creates a margin of 'margin_px' pixels on all four sides
ax = gca;
ax.Position = [margin_x, margin_y, 1 - 2*margin_x, 1 - 2*margin_y];

colors = lines(5);
markers = {'*', 's', 'd', '^', 'v'};

for i = 1:length(results_blockjacobi)
    hist = results_blockjacobi(i).history;
    flops = hist(:, 1);
    max_off = hist(:, 3);
    unique_sweeps = hist(:, 2);
    
    if i == 1
        loglog(flops, max_off, ...
        'LineWidth', 2, ...
        'LineStyle', ':', ...
        'Marker', markers{i}, ...
        'Color', colors(i,:), ...
        'DisplayName', results_blockjacobi(i).label);
    elseif i == 2
        loglog(flops, max_off, ...
        'LineWidth', 2, ...
        'LineStyle', '-.', ...
        'Marker', markers{i}, ...
        'Color', colors(i,:), ...
        'DisplayName', results_blockjacobi(i).label);
    else
        loglog(flops, max_off, ...
            'LineWidth', 2, ...
            'LineStyle', '-', ...
            'Marker', markers{i}, ...
            'Color', colors(i,:), ...
            'DisplayName', results_blockjacobi(i).label);
    end
    hold on;

    % Annotate final sweeps
    if i == 2
        text(1.1*flops(end), max_off(end), num2str(unique_sweeps(end)), ...
            'Color', colors(i,:), 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    elseif i == 4
        text(0.99*flops(end), 0.5 * max_off(end), num2str(unique_sweeps(end)), ...
        'Color', colors(i,:), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    elseif i == 1
        text(0.96*flops(end), 0.5 * max_off(end), num2str(unique_sweeps(end)), ...
        'Color', colors(i,:), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    elseif i == 3
        text(1.03*flops(end), 0.5 * max_off(end), num2str(unique_sweeps(end)), ...
        'Color', colors(i,:), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    else
        text(flops(end), 0.5 * max_off(end), num2str(unique_sweeps(end)), ...
        'Color', colors(i,:), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    end
end

% Decorations
yline(eps_threshold, '--r', ...
            'LineWidth', 1.5, ...
            'DisplayName', '1e-7 × max(max(abs(A)))');

xline(flops_limit, '--', ...
            'Color', [0.5, 0.5, 0.5], ...
            'LineWidth', 1.5, ...
            'DisplayName', 'flops limit');

xlabel('FLOPs');
ylabel('Max Abs Off-Diagonal');

grid on;
legend('Location', 'southwest'); legend('boxoff'); legend('AutoUpdate','off');
    set(legend, ...
    'EdgeColor', 'black', ...   % edge color
    'Color', 'white', ...       % background color
    'Box', 'on');               % display the edge

title('Block Jacobi Variants: flops vs. max off-diagonal');
