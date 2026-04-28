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

%% Q1a. Determine q[n] and frequency response Q(w)
% p[n] = h[n] (we don't care about T)
% q[n] = (p[n] * p*[-n]) / ||p||
%      = (h[n] * h*[-n]) / ||h||

mod_h = sqrt(sum(h_n.^2));
q_n = conv(h_n, fliplr(h_n)) / mod_h;
disp('q[n]:');
disp(q_n);

figure;
freqz(q_n, 1, 1024, 'whole');
title('Frequency Response |Q(w)|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');

%% Q1b. Best spectral efficiency for Pe = 1e-4, sigma2 = 0.1995
sigma2 = 0.1995
% We use the MFB for SNR here, since we're calculating our best case scenario. TODO: Check with Faraz or the Toddlers if this is right.
SNR_mfb_linear = (Es_bar * mod_h^2)/ sigma2; % ||p||^2 = ||h||^2
Gap_dB = 6.6; % From course reader
Gap_linear = 10^(Gap_dB/10);
b = log2(1 + SNR_mfb_linear / Gap_linear);

disp('Ideal spectral efficiency:')
disp(b)


