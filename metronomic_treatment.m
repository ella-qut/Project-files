% Metronomic treatment: Alternating periods of treatment and rest
function dydt = metronomic_treatment(t, y, r1, r2, w12, w21, delta)
    N1 = y(1);
    N2 = y(2);
    dummy = t;
    
    % Alternating between treatment and rest every 100 time steps
    if mod(t, 200) < 100
        % Treatment phase
        dydt(1) = (r1 - delta) * N1 - w12 * N1 + w21 * N2; % Growth of cancer type 1 with reduced rate
        dydt(2) = r2 * N2 - w21 * N2 + w12 * N1; % Growth of cancer type 2 with reduced rate
    else
        % Rest phase (no treatment)
        dydt(1) = r1 * N1 - w12 * N1 + w21 * N2;
        dydt(2) = r2 * N2 - w21 * N2 + w12 * N1;
    end
    
    dydt(3) = 0; % No change in PSA for simplicity
end