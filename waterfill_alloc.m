function p = waterfill_alloc(g, Pt)
    %% Inputs:
    % g  -> 1 x K vector of channel SNRs
    % Pt -> total power constraint
    
    % Outputs:
    % p  -> 1 x K vector of power allocations
    
    %% 0. Initialize variables
    gammas = 1 ./ g;              % Effective noise variances
    K = length(g);                % Total number of channels
    p = zeros(1, K);              % Final power allocation array
    active = true(1, K);          % Logical mask to track which channels are still on
    
    %% Waterfilling loop
    while true
        % 1. Calculate water level using active channels
        K_active = sum(active);
        water_level = (Pt + sum(gammas(active))) / K_active;
        
        % 2. Calculate power for active channels 
        p_temp = water_level - gammas(active);
        
        % 3. Check if any channels need to be turned off
        if any(p_temp < 0)
            % Identify the channels that failed in this round
            failed_in_temp = (p_temp < 0);
            
            % Map failuresto the 'active' mask to turn them off
            active_indices = find(active); 
            active(active_indices(failed_in_temp)) = false;
        else
            % If no powers are negative, waterfilling is complete.
            p(active) = p_temp; % Assign the positive powers to the final array
            break;
        end
    end
end


