%% Q2a.
clear all; close all; clc;
omega = linspace(0, pi, 1000);

q_amp = abs(sqrt(1.89)./(1.99 - 0.9*exp(1j*omega) + 0.8*exp(2j*omega) - 0.9*exp(-1j*omega) + 0.8*exp(-2j*omega)));


figure;
plot(omega, q_amp, 'LineWidth', 2);
grid on;
title('Amplitude of the Frequency Spectrum');
xlabel('\Omega (radians/sample)');
ylabel('|Q(e^{j\Omega})|');
xlim([0, pi]);
xticks([0, pi/2, pi]);
xticklabels({'0', '\pi/2', '\pi'});
%% Q2b.
Es = 1;

U = [0.58, -0.65, 1.37,  -0.65, 0.58,  0,     0,     0,     0;
     0,    0.58,  -0.65, 1.37,  -0.65, 0.58,  0,     0,     0;
     0,    0,     0.58,  -0.65, 1.37,  -0.65, 0.58,  0,     0;
     0,    0,     0,     0.58,  -0.65, 1.37,  -0.65, 0.58,  0;
     0,    0,     0,     0,     0.58,  -0.65, 1.37,  -0.65, 0.58];

Rw = (0.1/1.89) .* [1.89, -0.9, 0.8,  0,    0;
                    -0.9, 1.89, -0.9, 0.8,  0;
                    0.8,  -0.9, 1.89, -0.9, 0.8,
                    0,    0.8,  -0.9, 1.89, -0.9;
                    0,    0,    0.8,  -0.9, 1.89];
cmmse = inv(U * ctranspose(U) + (1/Es) * Rw) * U(:, 5)

%% Q2c.



