%% Q3
clear all; close all; clc;
load snr_vals.mat

% 1. Initialization
PtVecdB = -10:0.1:35; % SNR in dB
PtVecLin = 10.^(PtVecdB/10); % Convert SNR to linear scale
CapWfVec = zeros(1, length(PtVecLin)); % Vector of capacities when waterfilling is used
CapEqVec = zeros(1, length(PtVecLin)); % Vector of capacities when equal power allocation is used
CapHiVec = zeros(1, length(PtVecLin)); % Vector of capacities when the highest SNR channel is used

% 2. Calculate capacities for each SNR
for i = 1:length(PtVecLin)
    Pt = PtVecLin(i);
    p = waterfill_alloc(g, Pt);
    % 2a. Find capacity of waterfilled channel
    sum = 0;
    for j = 1:10
        sum = sum + log2(1 + p(j) * g(j));
    end
    CapWfVec(i) = sum;
    
    % 2b. Find capacity of equal power allocated channel
    sum = 0;
    for j = 1:10
        sum = sum + log2(1 + (Pt/10) * g(j));
    end
    CapEqVec(i) = sum;

    % 2c. Find capacity when only using max SNR channel
    CapHiVec(i) = log2(1 + Pt * max(g));
end

% 3. Plot everything
figure; 
plot(PtVecdB, CapWfVec)
hold on;
plot(PtVecdB, CapEqVec)
plot(PtVecdB, CapHiVec)
xlabel('(Pt) [dB]');
ylabel('Capacity [bits]');
title('Capacity vs Total Power');
legend('Water-filling', 'Equal Power', 'Max SNR Channel');

%% Q4
clear all; close all; clc;
T = 1; t = -10:0.01:10;

q = (1/sqrt(2.06)) * ( 0.5*sinc((t - 3*T) / T) - 0.45*sinc((t - 2*T) / T) - 0.9*sinc((t - T) / T) + 2.06*sinc(t / T) ...
                     + 0.5*sinc((t + 3*T) / T) - 0.45*sinc((t + 2*T) / T) - 0.9*sinc((t + T) / T));

% 0. Plot q
figure;
plot(t, q, 'b-', 'LineWidth', 1.5);
hold on;
intervals = min(t):T:max(t); 
xline(intervals, '--or', 'Alpha', 0.4);
grid on;
xlabel('Time (t)');
ylabel('Amplitude (q)');
title('Composite Signal with T-Interval Markers');

% 1. Define the sampling times
t_intersections = -10:T:10; 

% 2. Calculate 'q' at sample points
q = (1/sqrt(2.06)) * ( 0.5*sinc((t_intersections - 3*T) / T) - 0.45*sinc((t_intersections - 2*T) / T) - 0.9*sinc((t_intersections - T) / T) + 2.06*sinc(t_intersections / T) ...
                     + 0.5*sinc((t_intersections + 3*T) / T) - 0.45*sinc((t_intersections + 2*T) / T) - 0.9*sinc((t_intersections + T) / T));

% 3. Print values
for i = 1:length(t_intersections)
    if t_intersections(i) >= -3 && t_intersections(i) <= 3
        fprintf('%5d      |   %7.4f\n', t_intersections(i), q(i));
    end
end