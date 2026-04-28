function [E_n, b_n, b_bar] = ofdm_waterfilling(p, N, Gap, Es_bar, sigma2)
    % 1. Calculate parameters
    L = length(p) - 1;
    nu = L; % Optimum CP length
    E_tot = N * Es_bar; % Total transmit energy budget for OFDM
    
    % 2. Calculate channel gains via N-point FFT
    H = fft(p, N);
    g = (abs(H).^2) / sigma2;
    
    % 3. For SVD or OFDM, we divide our channel gains by the Gap
    effective_gammas = g / Gap;
    
    % 4. Call your pre-existing waterfilling function
    E_n = waterfill_alloc(effective_gammas, E_tot);
    
    % 5. Calculate rates and spectral efficiency
    b_n = zeros(1, N);
    active = E_n > 0;
    b_n(active) = 0.5 * log2(1 + (E_n(active) .* g(active)) / Gap);
    
    b_bar = sum(b_n) / (N + nu); % Spectral Efficiency
end