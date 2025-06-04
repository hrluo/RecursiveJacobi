function A = generate_spike_spectrum_matrix(n, num_spikes, spike_ratio)
% ------------------------------------------------------------------------
% Generates a symmetric matrix with spike spectra, meaning a few large
% eigenvalues (spikes) and many small ones. This results in a matrix with a
% large condition number.
%
% Inputs:
%   - n: Size of the matrix (n x n).
%   - num_spikes: Number of large eigenvalues (spikes).
%   - spike_ratio: Ratio of spike eigenvalues to small ones (spike/small).
%
% Output:
%   - A: Generated symmetric matrix with spike spectra.
% ------------------------------------------------------------------------

    % Step 1: Generate the eigenvalue spectrum
    small_eigenvalue = 1;  % Set a small eigenvalue for the majority
    large_eigenvalue = spike_ratio;  % Large eigenvalues (spikes)

    % Create a diagonal matrix of eigenvalues: few large, the rest small
    eigenvalues = small_eigenvalue * ones(n, 1);
    eigenvalues(1:num_spikes) = large_eigenvalue;  % Assign large eigenvalues (spikes)

    % Step 2: Generate a random orthogonal matrix (using QR decomposition)
    [Q, ~] = qr(randn(n, n));

    % Step 3: Create the symmetric matrix A with the desired spectrum
    A = Q * diag(eigenvalues) * Q';  % A = Q * Lambda * Q'
    
    % Ensure symmetry due to numerical errors
    A = (A + A') / 2;
end
