%% Q1f.
clear all; close all; clc;

% Parameters
alphaValues = [0.25, 0.5, 1]; % Excess bandwidth values
m = 32;                       % Oversampling factor (samples per symbol period)
span = 5;                     % Truncate at [-5T, 5T]

figure;
grid on;
hold on;
for alpha = alphaValues
    [waveform, t_norm] = raised_cosine(alpha, m, span); 
    
    % Plot waveform vs normalized time
    plot(t_norm, waveform, 'LineWidth', 1.5); 
end
hold off;

xlabel('Normalized Time (t/T)');
ylabel('Amplitude');
title('RCPs with Different Excess Bandwidths');
legend(arrayfun(@(x) sprintf('\\alpha = %.2f', x), alphaValues, 'UniformOutput', false));
axis([-5 5 -0.3 1.1]);

% We note that increasing alpha / excess bandwidth causes the sidelobes to attenuate faster.


%% Q1g and Q1h.
clear all; close all; clc;

% Parameters
T = 1; 
alpha = 0.5;
m = 16;               % Oversampling factor (samples per symbol period)
ts = T / m;           % Sampling interval

% FFT Setup
fs_desired = 1 / (64 * T);               % Desired frequency granularity
Nmin = ceil(1 / (fs_desired * ts));      % Minimum length DFT
Nfft = 2^(nextpow2(Nmin));               % FFT size (power of 2)

fs = 1 / (Nfft * ts);                    % Actual frequency resolution
% Normalized frequency vector (f * T)
freqs_norm = ((1:Nfft) - 1 - Nfft/2) * fs * T; 

% Generate time domain signals with the given RCP function
[rc_2, ~] = raised_cosine(alpha, m, 2); % Truncation at [-2T, 2T]
[rc_5, ~] = raised_cosine(alpha, m, 5); % Truncation at [-5T, 5T]

% Frequency domain signals
% fft function automatically zero-pads the time-domain signal to Nfft
RC_2_f = ts * fft(rc_2, Nfft);
RC_5_f = ts * fft(rc_5, Nfft);

% Center the spectrums
RC_2_f_centered = fftshift(RC_2_f);
RC_5_f_centered = fftshift(RC_5_f);

figure;
plot(freqs_norm, abs(RC_2_f_centered), 'r', 'LineWidth', 1.5);
hold on;
plot(freqs_norm, abs(RC_5_f_centered), 'b', 'LineWidth', 1.5);

xlabel('Normalized Frequency (fT)');
ylabel('Magnitude Spectrum');
title('Effect of Truncation on Spectrum');
legend('Truncation [-2T, 2T]', 'Truncation [-5T, 5T]', 'Location', 'best');
xlim([-2 2]);
grid on;
hold off;


Nfft_high = 2^20; 

% Recompute FFTs with zero-padding
RC_2_f_high = ts * fft(rc_2, Nfft_high);
RC_5_f_high = ts * fft(rc_5, Nfft_high);

% New high-res frequency vector
fs_actual = 1 / ts;
freqs_high_norm = ((1:Nfft_high) - 1 - Nfft_high/2) * fs_actual * T;

P2 = abs(fftshift(RC_2_f_high)).^2;
P5 = abs(fftshift(RC_5_f_high)).^2;

% sort the frequencies by their absolute distance from the center
[f_abs_sorted, sort_idx] = sort(abs(freqs_high_norm));
P2_sorted = P2(sort_idx);
P5_sorted = P5(sort_idx);

% Cumulative energy. Normalize to total energy (0 to 1)
cum_E2 = cumsum(P2_sorted);
cum_E5 = cumsum(P5_sorted);
norm_E2 = cum_E2 / cum_E2(end);
norm_E5 = cum_E5 / cum_E5(end);

% Calculate 95% energy bandwidth
target_95 = 0.95;
idx_95_2 = find(norm_E2 >= target_95, 1);
idx_95_5 = find(norm_E5 >= target_95, 1);
bw_95_2 = f_abs_sorted(idx_95_2);
bw_95_5 = f_abs_sorted(idx_95_5);

% Calculate 10^-8 energy bandwidth
target_high = 1 - 1e-8;
idx_high_2 = find(norm_E2 >= target_high, 1);
idx_high_5 = find(norm_E5 >= target_high, 1);
bw_high_2 = f_abs_sorted(idx_high_2);
bw_high_5 = f_abs_sorted(idx_high_5);

fprintf('Baseband Bandwidth Calculations (alpha = %.2f)\n', alpha);
fprintf('Nominal Ideal Bandwidth: %.3f / T\n', (1+alpha)/2);

fprintf('95%% Energy Containment Bandwidth:');
fprintf('Truncation [-2T, 2T]: %.4f / T\n', bw_95_2);
fprintf('Truncation [-5T, 5T]: %.4f / T\n', bw_95_5);
% We see that the pulse energy is closely matched in terms of containment
% bandwidth.

fprintf('(1 - 1e-8) Energy Containment Bandwidth:');
fprintf('Truncation [-2T, 2T]:      %.4f / T\n', bw_high_2);
fprintf('Truncation [-5T, 5T]:      %.4f / T\n', bw_high_5);
fprintf('(Note: -2T truncation forces you to look much further for the tails)\n');
% We see that truncating the pulse at +-2T causes the pulse energy to
% be spread much wider in frequency.