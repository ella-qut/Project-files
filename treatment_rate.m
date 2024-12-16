function [r, met_treatment_status] = treatment_rate(t, r1, r_treat, treat_time, no_treat_time)

cycle_length = treat_time + no_treat_time;
cycle_moment = mod(t, cycle_length);

if cycle_moment < treat_time
    r = r_treat;
    met_treatment_status = 1;

else
    r = r1;
    met_treatment_status = 0;
end

end