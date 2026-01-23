% Q5

clear; close all; clc;

N = 100;
Ns = 10; % Number of symbols to be plotted
Tsym = 1; % Symbol time in seconds
oversampleFactor = 500;
Tres = Tsym / oversampleFactor;
t = 0:Tres:Ns*Tsym-Tres; % Time axis for symbols to be plotted

bc = sign(2 * rand(1,N) - 1);
bs = sign(2 * rand(1,N) - 1);

uc = zeros(size(t));
us = zeros(size(t));

for n = 1:Ns
    idx = (t >= (n-1)*Tsym) & (t < n*Tsym);
    uc(idx) = bc(n);
    us(idx) = bs(n);
end

%% 5a. Plot baseband signals
figure
subplot(2,1,1)
plot(t, uc, 'LineWidth', 1.5)
ylim([-1.5 1.5])
grid on
title('u_c(t) over 10 symbols')
xlabel('t')
ylabel('Amplitude')

subplot(2,1,2)
plot(t, us, 'LineWidth', 1.5)
ylim([-1.5 1.5])
grid on
title('u_s(t) over 10 symbols')
xlabel('t')
ylabel('Amplitude')

%% 5b and 5c. Upconvert baseband real signal, 4 symbols
t = 0:Tres:4*Tsym-Tres;
Icarrier = cos(40*pi*t);
Qcarrier = sin(40*pi*t);

uc = zeros(size(t));
us = zeros(size(t));

for n = 1:4
    idx = (t >= (n-1)*Tsym) & (t < n*Tsym);
    uc(idx) = bc(n);
    us(idx) = bs(n);
end

up1 = uc .* Icarrier;
up2 = us .* Qcarrier;

up = up1 - up2;

figure
subplot(2,1,1)
plot(t, up1, 'LineWidth', 1.5)
grid on
title('u_{p,1}(t) over 4 symbols')
xlabel('t')
ylabel('Amplitude')

subplot(2,1,2)
plot(t, up2, 'LineWidth', 1.5)
grid on
title('u_{p,2}(t) over 4 symbols')
xlabel('t')
ylabel('Amplitude')

figure 
plot(t, up, 'LineWidth', 1.5)
grid on
title('u_p(t) over 4 symbols')
xlabel('t')
ylabel('Amplitude')

%% 5d. Demod, 10 symbols
t = 0:Tres:Ns*Tsym-Tres;
Icarrier = cos(40*pi*t);
Qcarrier = sin(40*pi*t);

uc = zeros(size(t));
us = zeros(size(t));

for n = 1:Ns
    idx = (t >= (n-1)*Tsym) & (t < n*Tsym);
    uc(idx) = bc(n);
    us(idx) = bs(n);
end

up1 = uc .* Icarrier;
up2 = us .* Qcarrier;

up = up1 - up2;


mixI = 2*up .* cos(40*pi*t);
mixQ = 2*up .* sin(40*pi*t);

h  = ones(1, 0.25 / Tres);

vc = conv(mixI, h, 'same') * Tres;
vs = conv(mixQ, h, 'same') * Tres;

figure
subplot(2,1,1)
plot(t, uc, 'LineWidth', 1.2); hold on
plot(t, vc, 'LineWidth', 1.2); grid on
title('I: u_c(t) vs v_c(t)'); xlabel('t'); ylabel('Amplitude')
legend('u_c(t)','v_c(t)')

subplot(2,1,2)
plot(t, us, 'LineWidth', 1.2); hold on
plot(t, vs, 'LineWidth', 1.2); grid on
title('Q: u_s(t) vs v_s(t)'); xlabel('t'); ylabel('Amplitude')
legend('u_s(t)','v_s(t)')