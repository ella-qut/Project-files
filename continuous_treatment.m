% Continuous treatment: Growth rate reduction
function dydt = continuous_treatment(t, y, r1, r2, w12, w21, delta)
    N1 = y(1);
    N2 = y(2);
    
    dydt = zeros(3, 1); % Pre-allocate the derivative vector
    
    dyN1 = (r1 - delta) * N1 - w12 * N1 + w21 * N2; % Growth of cancer type 1 with reduced rate
    dyN2 = r2 * N2 - w21 * N2 + w12 * N1; % Growth of cancer type 2 with reduced rate
    dyPSA = 0; % No change in PSA for simplicity
    dummy = t;
    dydt = [dyN1; dyN2; dyPSA];
end
