clear all; close all; clc;
addpath 'C:\Users\Ethan\OneDrive\Documents\MATLAB\DigiComm';
%% Q1a. Find singular values;
% Parameters
N = 8;
L = 1;
Esbar = 1;
Gap = 2;
T = 1; % Note that T does NOT affect singular values!


c = [1, 1, zeros(1, N-1)];
r = [1, zeros(1, N-1)];
P = 1/sqrt(T) * toeplitz(c,r);
lambdas = svd(P) % Lambdas correspond to channel gains, i.e. |h(n)|

%% Q1b. Waterfilling to find optimal powers and bitrates
sigma2 = 0.2;
pTotal = (N + L) * Esbar;
Gammas = ((lambdas .^ 2) / (sigma2))';
En = waterfill_alloc(Gammas / Gap, pTotal);
SNRn = En .* Gammas;
rn = 0.5 * log2(1 + SNRn ./ Gap);

svd_results = table((1:N)', En', rn', 'VariableNames', {'n', 'Power', 'Bit Rate'})


%% Q1c. 
% Find per subchannel rate
rBar = sum(rn) / (N + L);
SNRsvd = Gap * (2^(2 * rBar) - 1);
SNRsvd_dB = 10 * log10(SNRsvd);
fprintf('Multicarrier SNR for SVD: %0.2f dB \n', SNRsvd_dB)
% SNRmfb is 10dB. We are 1.32dB away from the matched filter bound!

%% Q3b: 1 + 0.05D Channel
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('3b: 1 + 0.05D Channel');
p1 = [1, 0.05];
N = 8;
Gap = 1;
Es_bar = 1;
sigma2 = 0.125;

[E_n, b_n, b_bar] = ofdm_waterfilling(p1, N, Gap, Es_bar, sigma2);

disp('Energy Distribution per subchannel (E_n):'); 
disp(E_n);
disp('Rates per subchannel (b_n):'); 
disp(b_n);
disp(['Spectral Efficiency (b_bar): ', num2str(b_bar), ' bits/dim']);
disp(' ');


%% Q3c: 1 + 0.5D + D^2 - D^3
disp('3c: 1 + 0.5D + D^2 - D^3 Channel');
p2 = [1, 0.5, 1, -1];
Gap = 2;
Es_bar = 1;
sigma2 = 0.1;


% Iterate over N = 2^k for k = 2 to 8
k_vals = 2:8;
N = 2.^k_vals;
b_bar_vals = zeros(size(N));

for i = 1:length(N)
    [~, ~, b_bar_vals(i)] = ofdm_waterfilling(p2, N(i), Gap, Es_bar, sigma2);
end

% Plot results!
figure;
plot(N, b_bar_vals, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Number of Subchannels (N)');
ylabel('Spectral Efficiency (bits/dim)');
title('OFDM Spectral Efficiency vs. N for 1 + 0.5D + D^2 - D^3');
grid on;

N_max = 256;
b_bar_max = b_bar_vals(end); % Last calculated value is for N=256

% Calculate Multicarrier SNR
% Formula: SNR_MC = Gamma * (2^(2*b_bar) - 1)
SNR_MC_linear = Gap * (2^(2 * b_bar_max) - 1);
SNR_MC_dB = 10 * log10(SNR_MC_linear);

% Calculate Matched Filter Bound (MFB) SNR
% Formula: SNR_MFB = (||p||^2 * Es_bar) / (N0/2)
norm_p2_sq = sum(p2.^2); % Sum of squared taps
SNR_MFB_linear = (norm_p2_sq * Es_bar) / sigma2;
SNR_MFB_dB = 10 * log10(SNR_MFB_linear);

disp(['Multicarrier SNR (N=256):     ', num2str(SNR_MC_dB), ' dB']);
disp(['Matched Filter Bound (MFB):   ', num2str(SNR_MFB_dB), ' dB']);
disp(['Gap from MFB:                 ', num2str(SNR_MFB_dB - SNR_MC_dB), ' dB']);





