%% FILES:
% main script: exp_updated.m
% RK4.m function file (ODE solver)
% treatment_rate.m function file (for metronomic treatment)
% adp_treat_rate.m function file (for adaptive treatment)

%% Establishing Variables/Growth Rates
r1 = 0.2; %Cell population 1 intrinsic growth rate
r2 = 0.2; %Cell population 2 intrinsic growth rate
w12 = 0.5; %Phenotypic switching rate from population 1 to 2 - potentially make switching rates higher (simulate parameter space)
w21 = 0.7; %Phenotypic switching rate from population 2 to 1
N1_0 = 100; %Initial population for population 1
N2_0 = 100; %Initial population for population 1

%Treatment model time boundaries (days):
a = 0; 
b = 600;
r_treat = -0.5; %Updated growth rate for sensitive population under treatment
n = 600; %number of time steps

% specific to metronomic model:
treat_time = 10; % time-step duration of treatment
no_treat_time = 30; % time-step duration of no treatment

% specific to adaptive model
PSA_0 = N1_0 + N2_0; %Initial PSA count
PSA_threshold = 0.5 * PSA_0; %Threshold for treatment stop/start

% graphing
y_max = 100;
log_y_max = log10(y_max);
%% Establishing ODEs for Cell population changes

% y = N1
% v = N2

% f1 = dN1/dt
% f2 = dN2/dt

f1 = @(t,y,v) r1*y - w12*y + w21*v;
f2 = @(t,y,v) r2*v - w21*v + w12*y;


%% Applying treatment methods to predict tumour behaviour:

treatment_method = 'adaptive'; % options: 'no_treatment', 'continuous', 'metronomic', 'adaptive'

if strcmp(treatment_method, 'no_treatment')
    f1 = @(t,y,v) r1*y - w12*y + w21*v;
    f2 = @(t,y,v) r2*v - w21*v + w12*y;
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    figure
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    hold on
    plot(t,v,'b', 'DisplayName', 'Resistant Population')
    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('No Treatment Exponential Model')
    ylim([-10,y_max])





elseif strcmp(treatment_method, 'continuous')
    f1 = @(t,y,v) r_treat*y - w12*y + w21*v;
    f2 = @(t,y,v) r2*v - w21*v + w12*y;
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    figure
    hold on
    set(gca, 'Color', [0.8,0.8,0.8]);
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    plot(t,v,'b', 'DisplayName', 'Resistant Population')

    treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
    set(treatment_box, 'DisplayName', 'Treatment Applied');

    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('Continuous Treatment Exponential Model')
    legend show;
    ylim([-10,y_max])






elseif strcmp(treatment_method, 'metronomic')
    treat_rate = @(t) treatment_rate(t, r1, r_treat, treat_time, no_treat_time); %function to determine treatment rate based on time in treatment cycle
    %e.g. whether treatment is being applied or not
    f1 = @(t,y,v) treat_rate(t)*y - w12*y + w21*v;
    f2 = @(t,y,v) r2*v - w21*v + w12*y;
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);

    met_treatment_status = false(size(t)); % Boolean array for treatment state

    for i = 1:length(t)
        [~, met_treatment_status(i)] = treatment_rate(t(i), r1, r_treat, treat_time, no_treat_time);
        %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

    end


    figure
    hold on
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    plot(t,v,'b', 'DisplayName', 'Resistant Population')

    x_limits = xlim;
    y_limits = ylim;

    met_treatment_regions = [0; met_treatment_status(:); 0];
    met_edges = diff(met_treatment_regions);
    met_start_idx = find(met_edges == 1);
    met_end_idx = find(met_edges == -1) -1;

    for i = 1:length(met_start_idx)
        met_x_start = t(met_start_idx(i));
        met_x_end = t(met_end_idx(i));        

        fill([met_x_start, met_x_end, met_x_end, met_x_start], [0, 0, y_max, y_max], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    end

    treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
    set(treatment_box, 'DisplayName', 'Treatment Applied');

    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('Metronomic Treatment Exponential Model')
    legend show;
    ylim([-10,y_max])




elseif strcmp(treatment_method, 'adaptive') %I still need to account for PSA Decay (dPSA/dt = N1 + N2 - 0.5*PSA)
    treat_rate_adp = @(t,y,v) adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold); %function to determine treatment rate based on time in treatment cycle
    %e.g. whether treatment is being applied or not
    f1 = @(t,y,v) treat_rate_adp(t,y,v)*y - w12*y + w21*v;
    f2 = @(t,y,v) r2*v - w21*v + w12*y;
    %f2 = @(t,y,v) treat_rate(t)*v - w21*v + w12*y; %I was just experimenting here!!! change treat_rate(t) back to r2
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);

    treatment_status = false(size(t)); % Boolean array for treatment state

    for i = 1:length(t)
        [~, treatment_status(i)] = adp_treat_rate(t(i), r1, r_treat, y(i), v(i), PSA_threshold);
        %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

    end

    figure
    hold on
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    plot(t,v,'b', 'DisplayName', 'Resistant Population')

    %set(gca, 'YScale', 'log')

    x_limits = xlim;
    y_limits = ylim;

    treatment_regions = [0; treatment_status(:); 0];
    edges = diff(treatment_regions);
    start_idx = find(edges == 1);
    end_idx = find(edges == -1) -1;

    for i = 1:length(start_idx)
        x_start = t(start_idx(i));
        x_end = t(end_idx(i));        

        fill([x_start, x_end, x_end, x_start], [0, 0, y_max, y_max], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    end

    treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
    set(treatment_box, 'DisplayName', 'Treatment Applied');

    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('Adaptive Treatment Exponential Model');
    legend show;
    ylim([-10,y_max])
    

else
    error('Invalid treatment method specified');
end 


%% Parameter Space Model:

% Establishing range for w12 and w21 rates:

w_start = 0;
w_end = 5;
w_int = 500;

w12_vals = linspace(w_start, w_end, w_int);
w21_vals = linspace(w_start, w_end, w_int);


% Establishing Point for Progression
prog_pnt = 1000000;


% Initialising matrix for w12 and w21 combinations:
param_space = zeros(length(w12_vals), length(w21_vals));


% Loop to test adaptive model on each combination of w12 and w21:
for i = 1:length(w12_vals)
    for j = 1:length(w21_vals)
        %Assign w12 and w21
        w12 = w12_vals(i);
        w21 = w21_vals(j);

        %Compute RK4
        treat_rate_adp = @(t,y,v) adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold); %function to determine treatment rate based on time in treatment cycle
        f1 = @(t,y,v) treat_rate_adp(t,y,v)*y - w12*y + w21*v;
        f2 = @(t,y,v) r2*v - w21*v + w12*y;
        [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);


        %Creating a vector to compute the cell population at each time step
        %for the particular w12 w21 combination:
        current_pop = y+v;

        %Finding time taken to reach the progression point:
        time_to_prog = find(current_pop > prog_pnt, 1); %finds the first occurence of the population exceeding the progression point
        if isempty(time_to_prog)
            param_space(i,j) = 1e4; %assigning a large time index if the population did not exceed the progression point within the time domain
        else
            param_space(i,j) = t(time_to_prog); %assigning the time value using the index found.
        end
    end

end

%Creating colour map
figure;
imagesc(w12_vals, w21_vals, param_space);
colorbar;
xlabel('w_{12} values');
ylabel('w_{21} values');
title('Time to progression (days)')
set(gca, 'YDir', 'normal'); 




%% NOTES:

%When playing with the treatment rate nothing seemed to be happening to the
%treatment vs. treatment.

%When i changed the rate of the second population to determine whether that
%was just growing to fast to reach below the threshold, alternation in
%treatment did occur. I believe the nature of the exponential model pay
%prevent adaptive therapy from occuring properly and is not appropriate to
%model the protocol.

%Was experimenting in line 84. Would normally change treat_rate(t) back to
%r2 but was just checking that the adaptive mechanism worked in the model
%in the case the total tumour size decreased.

%Still need to account for the PSA decay? dPSA/dt = N1 + N2 - 0.5*PSA
%Might update the RK4 to solve for this DE as well.