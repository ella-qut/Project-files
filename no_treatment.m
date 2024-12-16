% No treatment: Exponential growth with no reduction

function dydt = no_treatment(t, y, r1, r2, w12, w21)
    N1 = y(1);
    N2 = y(2);
    dydt = zeros(2, 1);
    dyN1 = r1 * N1 - w12 * N1 + w21 * N2;
    dyN2 = r2 * N2 - w21 * N2 + w12 * N1;
    dummy = t;

    dydt = [dyN1; dyN2];
end