function [l, treatment_status, direction] = adp_tr_log(t, l1, l1_treat, y, v, PSA_threshold_low, PSA_threshold_high)


direction = 0;
PSA = max(y+v,0);

if PSA > PSA_threshold_high
    l = l1_treat;
    treatment_status = 1;
    direction = 0;
elseif PSA < PSA_threshold_low
    l = l1;
    treatment_status = 0;
    direction = 1;

else %in between the high and low threshold
    if direction == 0
        l = l1_treat;
        treatment_status = 1;
    else
        l = l1;
        treatment_status = 0;
    end
end

end