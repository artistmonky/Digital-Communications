function [powers] = waterfill_alloc(gammas, pTotal)
    % Inputs:
    % gammas - 1 x K vector of channel gains, defined by gamma = (|h|^2)/sigma2
    % pTotal - Total allowable transmit power
    % sigma2 - Noise variance

    active = true(1, length(gammas)); % Mask that indicates which channels are active
    while true
        K = sum(active); % Number of active subchannels
        a = (pTotal + sum(1 ./ gammas(active))) / K;
        powers(active) = a - (1 ./ gammas(active));
        if all(powers >= 0)
            break;
        else
            % Create a temporary array for searching
            temp_powers = powers;
            
            % Set inactive channels to Infinity so min() ignores them
            temp_powers(~active) = Inf; 
            
            % Find the lowest power channel and drop it (guaranteed to be be active and have negative power)
            [~, min_idx] = min(temp_powers);
            active(min_idx) = false;
            powers(min_idx) = 0;
        end
    end
end


