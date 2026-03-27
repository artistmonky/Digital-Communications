clear all; close all; clc;

sigma = 1;
A = sqrt(2);

% Q-function
Q = @(x) 0.5*erfc(x/sqrt(2));

% 1. Create p array (avoid exactly 1 to prevent log(0))
p_vals = linspace(0.5, 0.999, 10000);
B_vals = zeros(size(p_vals));
Pe_PAM = zeros(size(p_vals));
Pe_OOK = zeros(size(p_vals));

% 2. Loop through p values
for k = 1:length(p_vals)
    p = p_vals(k);
    
    B = sqrt(1 / (2*(1 - p)));
    B_vals(k) = B;
    
    gamma = (sigma^2/(2*A)) * log(p/(1-p));
    Pe = p     * Q((gamma + A)/sigma) + ...
         (1-p) * Q((A - gamma)/sigma);
    Pe_PAM(k) = Pe;    

    gamma_OOK = B + (sigma^2/(2*B)) * log(p/(1-p));
    
    Pe_o = p * Q(gamma_OOK/sigma) + ...
           (1-p) * Q((2*B - gamma_OOK)/sigma);
    Pe_OOK(k) = Pe_o;
end

% 3. Graph probability of error
figure;
plot(p_vals, Pe_PAM, 'LineWidth', 2);
hold on;
plot(p_vals, Pe_OOK, 'LineWidth', 2);
grid on;
xlabel('p');
ylabel('Probability of Error');
title('Probability of Error vs p');

legend('2-PAM','OOK');

% 4. Find intersection
diffVals = Pe_PAM - Pe_OOK;
intersectionIndex = find(diffVals(1:end-1).*diffVals(2:end) < 0, 1);
p_intersection = p_vals(intersectionIndex)


