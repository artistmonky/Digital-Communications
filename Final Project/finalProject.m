clear all; close all; clc

% 0.1 Declare Channel Parameters
delay_ns = [0, 30, 70, 90, 110, 190,   410];
power_dB = [0, -1, -2, -3, -8,  -17.2, -20.8];
Fs = 10e6; % 10 MHz 
Ts = 1/Fs; % 100 ns
Es_bar = 1;


% 0.2 Find channel impulse response h[n]
tap_indices = round((delay_ns * 1e-9) / Ts);

% Initialize discrete-time channel pulse response
L = max(tap_indices); % Maximum delay index (check that L = 4)
h_n = zeros(1, L + 1); % We initialize h[n] with zeros and only fill in the taps with responses

% Assign tap values (later assignments overwrite earlier ones for same index)
for i = 1:length(delay_ns)
    idx = tap_indices(i) + 1; % 1-based indexing for MATLAB
    h_n(idx) = sqrt(10^(power_dB(i)/10)); % Linear voltage gain. 
end
disp('h[n]:');
disp(h_n); 
% We should see these taps (in linear units):
% n = 0 -> 0.891
% n = 1 -> 0.398
% n = 2 -> 0.138
% n = 3 -> 0
% n = 4 -> 0.0912
% TODO: Check with Toddlers if these values are kosher 
disp('---------------------------------------------------------------------')
disp(' ')

%% Q1a. Determine q[n] and frequency response Q(w)
% p[n] = h[n] (we don't care about T)
% q[n] = (p[n] * p*[-n]) / ||p||
%      = (h[n] * h*[-n]) / ||h||

mod_h = sqrt(sum(h_n.^2));
mod_p = mod_h;
q_n = conv(h_n, fliplr(h_n)) / mod_h;
disp('Q1a:')
disp('q[n]:');
disp(q_n);

figure;
freqz(q_n, 1, 1024);
title('Frequency Response |Q(w)|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
disp('---------------------------------------------------------------------')
disp(' ')

%% Q1b. Best spectral efficiency for Pe = 1e-4, sigma2 = 0.1995
sigma2 = 0.1995;
% We use the MFB for SNR here, since we're calculating our best case scenario. TODO: Check with Faraz or Toddlers if this is right.
SNR_mfb_linear = (Es_bar * mod_p^2)/ sigma2; % ||p||^2 = ||h||^2
Gap_dB = 6.6; % From course reader
Gap_linear = 10^(Gap_dB/10);
b = log2(1 + SNR_mfb_linear / Gap_linear);

disp('Q1b:')
disp('Ideal spectral efficiency:')
disp(b)
disp('---------------------------------------------------------------------')
disp(' ')

%% Q1c. Pe for 2PAM with SNR_dB = -5:2:11
% TODO: "we transmit 2PAM symbols over this channel without any
% equalization" - does this mean with or without MF? Check with Faraz
SNR_dB_vec = -5:2:11;
SNR_linear_vec = 10.^(SNR_dB_vec / 10);
I = Es_bar * (sum(q_n .^ 2) - (q_n(ceil(end/2)) ^ 2)); % TODO: If we're not supposed to use a MF, swap this for I = sum(h_n(2:end).^2) * Es_bar
S = mod_p^2 * Es_bar; % TODO: If we're not supposed to use a MF, swap mod_p^2 with h_n(1)^2
sigma2_vec = S ./ SNR_linear_vec; 

SINR_linear_vec = zeros(size(SNR_dB_vec));
for i = 1:length(SNR_dB_vec)
    SINR_linear_vec(i) = S / (I + sigma2_vec(i)); % Find SINR for each noise level
end

Pe_vec_estimated = zeros(size(SNR_dB_vec)); 
for i = 1:length(SNR_dB_vec)
    Pe_vec_estimated(i) = qfunc(sqrt(SINR_linear_vec(i)));
end

disp('Q1c:')
Pe_results = table(SNR_dB_vec', Pe_vec_estimated', 'VariableNames', {'SNR (dB)', 'Pe (using NNUB)'})
disp('---------------------------------------------------------------------')
disp(' ')

%% Q1d. Pe Analysis via Monte Carlo Simulation
num_iterations = 10000;
num_symbols = 1000;

% Q1d.1. No matched filtering
Pe_vec_empirical_noMF = zeros(size(SNR_dB_vec)); 
for i = 1:length(sigma2_vec) % 1. Loop through each noise level
    % Initialize error probability array for current noise level
    Pe_simulation = zeros(num_iterations, 1);
    
    for j = 1:num_iterations % 2. Loop through each iteration
        symbols = randi([0, 1], num_symbols, 1) * 2 - 1; % Generate random symbols for 2PAM. Map 0 to -1 and 1 to 1
        % Transmit symbols through the channel
        received_symbols = conv(symbols, h_n) + sqrt(sigma2_vec(i)) * randn(num_symbols + L, 1);
        received_symbols = received_symbols(1:num_symbols); % We have to get rid of the last L taps, as it's the just channel response from the last few symbols
        
        % Decode symbols
        detected_symbols = sign(received_symbols); % Threshold for 2-PAM is at 0
        
        % Calculate the number of errors
        Pe_simulation(j) = sum(detected_symbols ~= symbols);
    end
    
    % 3. Average Pe across iterations
    Pe_avg = mean(Pe_simulation) / num_symbols;
    % Store the average Pe for the current noise level
    Pe_vec_empirical_noMF(i) = Pe_avg; 
end


disp('Q1d:')
Pe_results = table(SNR_dB_vec', Pe_vec_empirical_noMF', 'VariableNames', {'SNR (dB)', 'Pe (empirical, no MF)'})
disp('---------------------------------------------------------------------')
disp(' ') 

% Q1d.2. With Matched filtering
% Recall that p[n] = h[n]. Thus:
% g_MF = p*[-n] / ||p||
%      = h*[-n] / ||h||
g_MF = fliplr(h_n) / mod_h; % NOTE! THIS FLIPS IT BUT THEN PUSHES THE FIRST INDEX TO t = 0. This MF is CAUSAL. We'll need to sample starting from L + 1.

Pe_vec_empirical_MF = zeros(size(SNR_dB_vec)); 
for i = 1:length(sigma2_vec) % 1. Loop through each noise level
    % Initialize error probability array for current noise level
    Pe_simulation_MF = zeros(num_iterations, 1);
    
    for j = 1:num_iterations % 2. Loop through each iteration
        symbols = randi([0, 1], num_symbols, 1) * 2 - 1; % Generate random symbols for 2PAM. Map 0 to -1 and 1 to 1
        
        received_symbols_full = conv(symbols, h_n) + sqrt(sigma2_vec(i)) * randn(num_symbols + L, 1);
        
        % Pass through Matched Filter
        mf_output = conv(received_symbols_full, g_MF);
        
        %Start sampling at L + 1.
        sampled_mf_output = mf_output(L + 1 : L + num_symbols);
        
        detected_symbols_MF = sign(sampled_mf_output);
        
        Pe_simulation_MF(j) = sum(detected_symbols_MF ~= symbols);     
    end
    
    % 3. Average Pe across iterations
    Pe_avg = mean(Pe_simulation_MF) / num_symbols;
    % Store the average Pe for the current noise level
    Pe_vec_empirical_MF(i) = Pe_avg; 
end


disp('Q1d:')
Pe_results = table(SNR_dB_vec', Pe_vec_empirical_MF', 'VariableNames', {'SNR (dB)', 'Pe (empirical, with MF)'})
disp('---------------------------------------------------------------------')
disp(' ') 

figure
plot(SNR_dB_vec, Pe_vec_estimated, 'k-', 'LineWidth', 2);
hold on;
plot(SNR_dB_vec, Pe_vec_empirical_noMF, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(SNR_dB_vec, Pe_vec_empirical_MF, 'rx--', 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;

grid on;
xlabel('SNR (dB)');
ylabel('Probability of Error (Pe)');
title('Comparison of Theoretical vs Empirical Error Probabilities');
legend('Theoretical (NNUB + MF)', 'Empirical (No MF)', 'Empirical (With MF)', 'Location', 'best');

% disp('Hi!! If this doesnt appear something is wrong')

% Q: Why there is a gap between the theoretical and simulated results? 
% A:

% Q: Under what conditions is the use of a matched filter beneficial? Provide an intuitive justification
% A: 

