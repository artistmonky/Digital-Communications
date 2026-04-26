clear all; close all; clc;

Es = 1;         
sigma2 = 0.1;   
T = 1;

% Q(D) = 1.89 - 0.9D - 0.9D^-1 + 0.8D^2 + 0.8D^-2
Q_k = [0.8, -0.9, 1.89, -0.9, 0.8];

%% Q2a.

w = linspace(0, pi, 1000);

q_amp = abs(sqrt(1.89)./(1.99 - 0.9*exp(1j*w) + 0.8*exp(2j*w) - 0.9*exp(-1j*w) + 0.8*exp(-2j*w)));

figure;
plot(w, q_amp, 'LineWidth', 2);
grid on;
title('Amplitude of the Frequency Spectrum');
xlabel('\Omega (radians/sample)');
ylabel('|Q(e^{j\Omega})|');
xlim([0, pi]);
xticks([0, pi/2, pi]);
xticklabels({'0', '\pi/2', '\pi'});

%% Q2b.
q_tilde = Q_k / sqrt(1.89);

% For a 5x9 matrix U using toeplitz:
% - The first column c defines the leftmost elements (padded with zeros)
% - The first row r defines the top elements (the full impulse response padded)
c = [q_tilde(1); 0; 0; 0; 0]; 
r = [q_tilde, 0, 0, 0, 0];
U = toeplitz(c, r);

% Noise covariance matrix
rw_col = (sigma2 / 1.89) * [1.89, -0.9, 0.8, 0, 0]'; 
Rww = toeplitz(rw_col);

% b[0] is at column 5 of U
u_0 = U(:, 5); 
c_mmse = (U * U' + (1/Es) * Rww) \ u_0

%% Q2c.

% Trruncate U to only include columns 3 to 7
U_trunc = U(:, 3:7);

% b[0] is the 3rd column of the truncated matrix
u_0_trunc = U_trunc(:, 3); 

c_mmse_trunc = (U_trunc * U_trunc' + (1/Es) * Rww) \ u_0_trunc

%% Q2d.

calc_sinr = @(c) ...
    (abs(c' * U(:, 5))^2) / ...                                 % Signal Power
    (sum(abs(c' * U(:, [1:4, 6:9])).^2) + real(c' * Rww * c));  % Interference + Noise Power

% Calculate SINR for both designs
sinr_b = calc_sinr(c_mmse);
sinr_c = calc_sinr(c_mmse_trunc);

fprintf('SINR for MMSE Equalizer (Part B):        %f \n', sinr_b);
fprintf('SINR for truncatedEqualizer (Part C):    %f\n', sinr_c);
% The SINR degrades when using the truncated equalizer.
% This is because the filter was optimized assuming b[-4], b[-3],
% b[3], and b[4] were zero, but they are present in the actual received
% signal and thus become unmitigated interference.

%% Q3c.
clear all;

% Parameters
SNR_mfb_dB = 10;
SNR_mfb = 10^(SNR_mfb_dB/10);
Es = 1;
Gamma = 2; 

% Noise variance calculation
sigma2 = (Es * 2) / SNR_mfb;

fprintf('N = 8\n');
N8 = 8;
P_total_8 = N8 * Es; % Total energy budget
n8 = 0:(N8-1);
g8 = 2 + 2*cos(2*pi*n8/N8); % Channel gains

[P8, lambda8] = waterfill_alloc(g8, P_total_8, Gamma, sigma2);

% Calculate Metrics
SNR_n_8 = (P8 .* g8) / sigma2;
r_n_8 = max(0, 0.5 * log2(1 + SNR_n_8 / Gamma));
r_bar_8 = sum(r_n_8) / N8;
SNR_mc_8 = Gamma * (2^(2 * r_bar_8) - 1);
SNR_mc_8_dB = 10 * log10(SNR_mc_8);
gap_to_mfb_8 = SNR_mfb_dB - SNR_mc_8_dB;

fprintf('n \t Gain (g_n) \t Power (P_n) \t SNR_n \t\t Rate (bits)\n');
fprintf('----------------------------------------------------------\n');
for i = 1:N8
    fprintf('%d \t %8.4f \t %8.4f \t %8.4f \t %8.4f\n', ...
        n8(i), g8(i), P8(i), SNR_n_8(i), r_n_8(i));
end
fprintf('\nMetrics:\n');
fprintf('Spectral Efficiency:         %.4f bits/symbol\n', r_bar_8);
fprintf('Multicarrier SNR:            %.4f dB\n', SNR_mc_8_dB);
fprintf('Distance to MFB:             %.4f dB\n\n', gap_to_mfb_8);


% Run for N = 16
fprintf('N = 16\n');
N16 = 16;
P_total_16 = N16 * Es;
n16 = 0:(N16-1);
g16 = 2 + 2*cos(2*pi*n16/N16);

% Perform Waterfilling
[P16, lambda16] = waterfill_alloc(g16, P_total_16, Gamma, sigma2);

% Calculate Metrics
SNR_n_16 = (P16 .* g16) / sigma2;
r_n_16 = max(0, 0.5 * log2(1 + SNR_n_16 / Gamma));
r_bar_16 = sum(r_n_16) / N16;
SNR_mc_16 = Gamma * (2^(2 * r_bar_16) - 1);
SNR_mc_16_dB = 10 * log10(SNR_mc_16);
gap_to_mfb_16 = SNR_mfb_dB - SNR_mc_16_dB;

% Print Tabular Results
fprintf('n \t Gain (g_n) \t Power (P_n) \t SNR_n \t\t Rate (bits)\n');
fprintf('----------------------------------------------------------\n');
for i = 1:N16
    fprintf('%d \t %8.4f \t %8.4f \t %8.4f \t %8.4f\n', ...
        n16(i), g16(i), P16(i), SNR_n_16(i), r_n_16(i));
end
fprintf('\nMetrics:\n');
fprintf('Spectral Efficiency:         %.4f bits/symbol/dim\n', r_bar_16);
fprintf('Multicarrier SNR:            %.4f dB\n', SNR_mc_16_dB);
fprintf('Distance to MFB:             %.4f dB\n\n', gap_to_mfb_16);

% Conclusion: Spectral Efficiency went UP from N=8 to N=16 with waterfilling


% Helper Function: Waterfilling Algorithm
function [P, lambda] = waterfill_alloc(g, P_total, Gamma, sigma2)
    N_ch = length(g);
    P = zeros(1, N_ch);
    
    % Only consider channels with strictly positive gain to avoid division by zero
    active_idx = find(g > 1e-6); 
    
    while true
        % Calculate the "inverse SNR" cost term for each active channel
        cost = (Gamma * sigma2) ./ g(active_idx);
        
        % Calculate the water level lambda
        lambda = (P_total + sum(cost)) / length(active_idx);
        
        % Tentative power allocation
        P_tentative = lambda - cost;
        
        % Check if any allocated power is negative
        if all(P_tentative >= 0)
            P(active_idx) = P_tentative;
            break; % Success! All active channels have non-negative power.
        else
            % Find the channel that was allocated the most negative power
            [~, min_idx] = min(P_tentative);
            dropped_ch = active_idx(min_idx);
            
            % Remove this channel from the active set and repeat
            active_idx(active_idx == dropped_ch) = [];
        end
    end
end



















% %% Q2a.
% clear all; close all; clc;
% omega = linspace(0, pi, 1000);
% 
% q_amp = abs(sqrt(1.89)./(1.99 - 0.9*exp(1j*omega) + 0.8*exp(2j*omega) - 0.9*exp(-1j*omega) + 0.8*exp(-2j*omega)));
% 
% 
% figure;
% plot(omega, q_amp, 'LineWidth', 2);
% grid on;
% title('Amplitude of the Frequency Spectrum');
% xlabel('\Omega (radians/sample)');
% ylabel('|Q(e^{j\Omega})|');
% xlim([0, pi]);
% xticks([0, pi/2, pi]);
% xticklabels({'0', '\pi/2', '\pi'});
% %% Q2b.
% 
% 
% 
% 
% Es = 1;
% sigma_sq = 0.1
% 
% U = [0.58, -0.65, 1.37,  -0.65, 0.58,  0,     0,     0,     0;
%      0,    0.58,  -0.65, 1.37,  -0.65, 0.58,  0,     0,     0;
%      0,    0,     0.58,  -0.65, 1.37,  -0.65, 0.58,  0,     0;
%      0,    0,     0,     0.58,  -0.65, 1.37,  -0.65, 0.58,  0;
%      0,    0,     0,     0,     0.58,  -0.65, 1.37,  -0.65, 0.58];
% 
% Rw = (0.1/1.89) .* [1.89, -0.9, 0.8,  0,    0;
%                     -0.9, 1.89, -0.9, 0.8,  0;
%                     0.8,  -0.9, 1.89, -0.9, 0.8,
%                     0,    0.8,  -0.9, 1.89, -0.9;
%                     0,    0,    0.8,  -0.9, 1.89];
% cmmse = inv(U * ctranspose(U) + (1/Es) * Rw) * U(:, 5)
% 
% %% Q2c. % DOUBLE CHECK WITH FARAZ THAT THIS IS WHAT YOU'RE SUPPOSED TO DO!
% % 1. Trim U to only keep b[-2] to b[2]
% U_trunc = [1.37,  -0.65, 0.58,  0,    0;     
%           -0.65, 1.37,  -0.65, 0.58,  0;
%           0.58,  -0.65, 1.37,  -0.65, 0.58;
%           0,     0.58,  -0.65, 1.37,  -0.65;
%           0,     0,     0.58,  -0.65, 1.37];
% 
% Rw = (0.1/1.89) .* [1.89, -0.9, 0.8,  0,    0;
%                     -0.9, 1.89, -0.9, 0.8,  0;
%                     0.8,  -0.9, 1.89, -0.9, 0.8,
%                     0,    0.8,  -0.9, 1.89, -0.9;
%                     0,    0,    0.8,  -0.9, 1.89];
% cmmse_trunc = inv(U_trunc * ctranspose(U_trunc) + (1/Es) * Rw) * U_trunc(:, 3)
% 
% %% Q2d. 
% 
% % Find SINR for original MMSE
% indices = 1:5;
% 
% % 1. Find Signal Energy
% u_0 = U(:, 5); 
% S = abs(cmmse' * u_0)^2;
% 
% % 2. Find Interference Energy
% U_interfere = U;
% U_interfere(:, 5) = 0; % just the interfering symbols
% I = sum(abs(cmmse' * U_interfere).^2);
% 
% % 3. Find Noise energy
% N = 0;
% for i = indices
%     ci = cmmse(i);
%     N = N + ci * sigma_sq;
% end
% 
% SINR = S / (I + N)




