%% Establishing Variables/Growth Rates
r1 = 0.15; %Cell population 1 intrinsic growth rate
r2 = 0.15; %Cell population 2 intrinsic growth rate
w12 = 0.3; %Phenotypic switching rate from population 1 to 2 - potentially make switching rates higher (simulate parameter space)
w21 = 0.6; %Phenotypic switching rate from population 2 to 1
N1_0 = 1000; %Initial population for population 1
N2_0 = 1000; %Initial population for population 1

%Treatment model time boundaries (days):
a = 0; 
b = 1000;
r_treat = -0.5; %Updated growth rate for sensitive population under treatment

n = 1000; %number of time steps
%n = 10000;

% specific to metronomic model:
treat_time = 10; % time-step duration of treatment
no_treat_time = 20; % time-step duration of no treatment

% specific to adaptive model
PSA_0 = N1_0 + N2_0; %Initial PSA count
PSA_threshold_low = 0.5 * PSA_0; %Threshold for treatment stop/start
PSA_threshold_high = 0.8* PSA_0;

% graphing
y_max = 1000000;
log_y_max = log10(y_max);

% Establishing Point for Progression
%prog_pnt = 1000000;
prog_pnt = 10000;

%% Applying treatment methods to predict tumour behaviour:
% No treatment
f1 = @(t,y,v) r1*y - w12*y + w21*v;
f2 = @(t,y,v) r2*v - w21*v + w12*y;
%[t_none,y_none,v_none] = RK4(f1,f2,a,b,n,N1_0,N2_0);
[t_none,y_none,v_none] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);


%%

% Continuous
f1 = @(t,y,v) r_treat*y - w12*y + w21*v;
f2 = @(t,y,v) r2*v - w21*v + w12*y;
%[t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
[t_con,y_con,v_con] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);


%%
% Metronomic
treat_rate = @(t) treatment_rate(t, r1, r_treat, treat_time, no_treat_time); %function to determine treatment rate based on time in treatment cycle
%e.g. whether treatment is being applied or not
f1 = @(t,y,v) treat_rate(t)*y - w12*y + w21*v;
f2 = @(t,y,v) r2*v - w21*v + w12*y;
%[t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
[t_met,y_met,v_met] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);

met_treatment_status = false(size(t_met)); % Boolean array for treatment state

for i = 1:length(t_met)
    [~, met_treatment_status(i)] = treatment_rate(t_met(i), r1, r_treat, treat_time, no_treat_time);
    %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

end


%%
% Adaptive

y_max = 10000;
treat_rate_adp = @(t,y,v) adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold_low, PSA_threshold_high); %function to determine treatment rate based on time in treatment cycle
%e.g. whether treatment is being applied or not
f1 = @(t,y,v) treat_rate_adp(t,y,v)*y - w12*y + w21*v;
f2 = @(t,y,v) r2*v - w21*v + w12*y;
%f2 = @(t,y,v) treat_rate(t)*v - w21*v + w12*y; %I was just experimenting here!!! change treat_rate(t) back to r2
%[t,y,v,stop_time] = RK4(f1,f2,a,b,n,N1_0,N2_0,prog_pnt);
%[t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
[t_adp,y_adp,v_adp] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);
total_count = y_adp + v_adp;
treatment_status = false(size(t_adp)); % Boolean array for treatment state

for i = 1:length(t_adp)
    [~, treatment_status(i), direction(i)] = adp_treat_rate(t_adp(i), r1, r_treat, y_adp(i), v_adp(i), PSA_threshold_low, PSA_threshold_high);
    %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

end

%% Time vs Pop Plots

tiledlayout(2,2);

% Tile 1
nexttile
plot(t_none,y_none,'r', 'DisplayName', 'Sensitive Population')
hold on
plot(t_none,v_none,'b', 'DisplayName', 'Resistant Population')
xlabel('Time, (days)')
ylabel('Cell Count')
title('No Treatment Exponential Model')
ylim([-10,y_max*100])
xlim([-1,200])



% Tile 2
nexttile
hold on
set(gca, 'Color', [0.8,0.8,0.8]);
plot(t_con,y_con,'r', 'DisplayName', 'Sensitive Population')
plot(t_con,v_con,'b', 'DisplayName', 'Resistant Population')

treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
set(treatment_box, 'DisplayName', 'Treatment Applied');

xlabel('Time, (days)')
ylabel('Cell Count')
title('Continuous Treatment Exponential Model')
legend show;
%ylim([-10,y_max])
xlim([-1,200])



% Tile 3
nexttile
hold on
plot(t_met,y_met,'r', 'DisplayName', 'Sensitive Population')
plot(t_met,v_met,'b', 'DisplayName', 'Resistant Population')

x_limits = xlim;
y_limits = ylim;

met_treatment_regions = [0; met_treatment_status(:); 0];
met_edges = diff(met_treatment_regions);
met_start_idx = find(met_edges == 1);
met_end_idx = find(met_edges == -1) -1;

for i = 1:length(met_start_idx)
    met_x_start = t_met(met_start_idx(i));
    met_x_end = t_met(met_end_idx(i));        

    fill([met_x_start, met_x_end, met_x_end, met_x_start], [0, 0, y_max, y_max], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
set(treatment_box, 'DisplayName', 'Treatment Applied');

xlabel('Time, (days)')
ylabel('Cell Count')
title('Metronomic Treatment Exponential Model')
legend show;
ylim([-10,10000])
xlim([-1,200])



% Tile 4
nexttile
hold on
plot(t_adp,y_adp,'r', 'DisplayName', 'Sensitive Population')
plot(t_adp,v_adp,'b', 'DisplayName', 'Resistant Population')
plot(t_adp,total_count, 'g', 'DisplayName', 'Total')

x_limits = xlim;
y_limits = ylim;

treatment_regions = [0; treatment_status(:); 0];
edges = diff(treatment_regions);
start_idx = find(edges == 1);
end_idx = find(edges == -1); %-1;

for i = 1:length(start_idx)
    x_start = t_adp(start_idx(i));
    % x_end = t(end_idx(i));

    % Determining x_end to avoid indexing errors:
    if end_idx(i) > length(t_adp)
        x_end = t_adp(end); % Assign the last value of t if index is out of bounds
    else
        x_end = t_adp(end_idx(i)); % Assign the valid index value
    end

    fill([x_start, x_end, x_end, x_start], [0, 0, y_max, y_max], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end

treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
set(treatment_box, 'DisplayName', 'Treatment Applied');

xlabel('Time, (days)')
ylabel('Cell Count')
title('Adaptive Treatment Exponential Model');
legend show;
ylim([-10,2000])
xlim([-1,200])
    



%% Parameter space plots:

tiledlayout(3,3);

r_start = 0.2;
r_end = 1;
r_int = 3;
r_vals = linspace(r_start, r_end, r_int);

rt_start = -0.1;
rt_end = -0.5;
rt_vals = linspace(rt_start, rt_end, r_int);

w_start = 0;
w_end = 1;
w_int = 100;

w12_vals = linspace(w_start, w_end, w_int);
w21_vals = linspace(w_start, w_end, w_int);


for g = 1:length(rt_vals)
    for h = 1:length(r_vals)

        r_treat = rt_vals(g);
        r2 = r_vals(h);
    
        % Initialising matrix for w12 and w21 combinations:
        param_space = zeros(length(w12_vals), length(w21_vals));
        
        
        % Loop to test adaptive model on each combination of w12 and w21:
        for i = 1:length(w12_vals)
            for j = 1:length(w21_vals)
                %Assign w12 and w21
                w12 = w12_vals(i);
                w21 = w21_vals(j);
        
                %Compute RK4
                treat_rate_adp = @(t,y,v) adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold_low, PSA_threshold_high); %function to determine treatment rate based on time in treatment cycle
                f1 = @(t,y,v) treat_rate_adp(t,y,v)*y - w12*y + w21*v;
                f2 = @(t,y,v) r2*v - w21*v + w12*y;
                [t,y,v] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);
        
        
                %Creating a vector to compute the cell population at each time step
                %for the particular w12 w21 combination:
                current_pop = y+v;
        
                %Finding time taken to reach the progression point:
                %time_to_prog = find(current_pop > prog_pnt, 1); %finds the first occurence of the population exceeding the progression point
                time_to_prog = find(movmean(current_pop,3)>prog_pnt, 1);
                if isempty(time_to_prog)
                    param_space(i,j) = 1000; %assigning a large time index if the population did not exceed the progression point within the time domain
                else
                    param_space(i,j) = t(time_to_prog); %assigning the time value using the index found.
                end
            end
        
        end
        
        %Creating colour map
        param_space_fixed = param_space';
        
        nexttile
        imagesc(w12_vals, w21_vals, param_space_fixed);
        colorbar;
        xlabel('w_{12} values');
        ylabel('w_{21} values');
        title(sprintf('r_{treat}=%.2f, r_2=%.2f', r_treat, r2))
        set(gca, 'YDir', 'normal');
    end
end