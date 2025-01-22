function [l, met_treatment_status] = met_tr_log(t, l1, l1_treat, treat_time, no_treat_time)

cycle_length = treat_time + no_treat_time;
cycle_moment = mod(t, cycle_length);

if cycle_moment < treat_time
    l = l1_treat;
    met_treatment_status = 1;

else
    l = l1;
    met_treatment_status = 0;
end

end