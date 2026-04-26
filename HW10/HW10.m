clear all; close all; clc;

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