function [Ts] = find_Ts(v, N)
    clc;
    Ts = (3 * N) / ((N + v) * 120e6)
    if (v * Ts) > 1e-6
        disp("Works")
    else
        disp("Fails")
    end
end