% Adaptive therapy: Adjust growth rate based on PSA levels
function dydt = adaptive_treatment(t, y, r1, r2, w12, w21, PSA_threshold_high, PSA_threshold_low, delta)
    N1 = y(1);
    N2 = y(2);
    PSA = y(3);
    dummy = t;
    dydt = zeros(3,1);

    % Determine treatment status based on PSA level
    if PSA > PSA_threshold_high
        % PSA is high, apply treatment (reduce growth rate)
        fprintf('Treatment ON at t = %.2f, PSA = %.2f\n', t, PSA);
        dyN1 = (r1 - delta) * N1 - w12 * N1 + w21 * N2; % Reduced growth during treatment
        dyN2 = r2 * N2 - w21 * N2 + w12 * N1;
    elseif PSA < PSA_threshold_low
        fprintf('Treatment OFF at t = %.2f, PSA = %.2f\n', t, PSA);
        % PSA is low, withdraw treatment (normal growth)
        dyN1 = r1 * N1; % Normal growth without treatment
        dyN2 = r2 * N2;
    elseif PSA < 0
        PSA = 0;
    else
        % Intermediate state: continue with reduced growth
        dyN1 = r1 * N1 - w12 * N1 + w21 * N2;
        dyN2 = r2 * N2 - w21 * N2 + w12 * N1;
    end
    
    % Update PSA level (assumes 1 cell = 1 PSA unit)
    dyPSA = (N1 + N2) - 0.5 * PSA; % PSA decay model

    % Constrain variables to avoid runaway growth
    N1 = min(N1, 1e8);
    N2 = min(N2, 1e8);
    PSA = min(PSA, 1e8);

    dydt = [dyN1; dyN2; dyPSA];

end
