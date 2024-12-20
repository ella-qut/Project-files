function [r, treatment_status] = adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold)
%PSA = y + v;

    PSA = max(y+v,0); % ensures PSA value does not go below 0
    %k = size(t);

    if PSA < PSA_threshold
        %disp('no treatment')
        r = r1;
        treatment_status = 0;
    else %PSA > PSA_threshold
        %disp('treatment')
        r = r_treat;
        treatment_status = 1;
    end

end