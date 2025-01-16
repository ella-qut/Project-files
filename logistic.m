%% Establishing Variables/growth rates

r1 = 1; %Cell population 1 intrinsic growth rate
r2 = 1; %Cell population 2 intrinsic growth rate
w12 = 0.02; %Phenotypic switching rate from population 1 to 2 - potentially make switching rates higher (simulate parameter space)
w21 = 0.03; %Phenotypic switching rate from population 2 to 1
l1 = 0.005; %Intensity of competition from population 2 onto population 1
l2 = 0.005; %Intensity of competition from population 1 onto population 2
N1_0 = 10; %Initial population for population 1
N2_0 = 10; %Initial population for population 1


y_max = 400;

%Treatment model time boundaries (days)
a = 0;
b = 400;
r_treat = -0.2;

n = 1800; %number of time steps

% specific to metronomic model:
treat_time = 100; % time-step duration of treatment
no_treat_time = 200; % time-step duration of no treatment

% specific to adaptive model
PSA_0 = N1_0 + N2_0; %Initial PSA count
PSA_threshold_low = 0.5 * PSA_0; %Threshold for treatment stop/start
PSA_threshold_high = 0.8* PSA_0;

prog_pnt = 300;

%% Applying treatment methods:

treatment_method = 'adaptive'; % options: 'no_treatment', 'continuous', 'metronomic', 'adaptive'

if strcmp(treatment_method, 'no_treatment')
    f1 = @(t,y,v) r1*y - w12*y + w21*v - l1*(y.^2 + y.*v);
    f2 = @(t,y,v) r2*v - w21*v + w12*y - l2*(v.^2 + y.*v);
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    %[t,y,v] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);
    figure
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    hold on
    plot(t,v,'b', 'DisplayName', 'Resistant Population')
    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('No Treatment Exponential Model')
    ylim([-10,y_max])





elseif strcmp(treatment_method, 'continuous')
    f1 = @(t,y,v) r1*y - w12*y + w21*v - l1*(y.^2 + y.*v);
    f2 = @(t,y,v) r2*v - w21*v + w12*y - l2*(v.^2 + y.*v);
    [t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    %[t,y,v] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);
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
    f1 = @(t,y,v) r1*y - w12*y + w21*v - l1*(y.^2 + y.*v);
    f2 = @(t,y,v) r2*v - w21*v + w12*y - l2*(v.^2 + y.*v);
    %[t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    [t,y,v] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);

    met_treatment_status = false(size(t)); % Boolean array for treatment state

    for i = 1:length(t)
        [~, met_treatment_status(i)] = treatment_rate(t(i), r1, r_treat, treat_time, no_treat_time);
        %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

    end


    figure
    hold on
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    plot(t,v,'b', 'DisplayName', 'Resistant Population')

    % semilogy(t,y,'r', 'DisplayName', 'Sensitive Population')
    % semilogy(t,v,'b', 'DisplayName', 'Resistant Population')

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
    treat_rate_adp = @(t,y,v) adp_treat_rate(t, r1, r_treat, y, v, PSA_threshold_low, PSA_threshold_high); %function to determine treatment rate based on time in treatment cycle
    %e.g. whether treatment is being applied or not
    f1 = @(t,y,v) r1*y - w12*y + w21*v - l1*(y.^2 + y.*v);
    f2 = @(t,y,v) r2*v - w21*v + w12*y - l2*(v.^2 + y.*v);
    %f2 = @(t,y,v) treat_rate(t)*v - w21*v + w12*y; %I was just experimenting here!!! change treat_rate(t) back to r2
    %[t,y,v,stop_time] = RK4(f1,f2,a,b,n,N1_0,N2_0,prog_pnt);
    %[t,y,v] = RK4(f1,f2,a,b,n,N1_0,N2_0);
    [t,y,v] = euler_stabCheck(f1,f2,a,b,n,N1_0,N2_0);
    total_count = y + v;
    treatment_status = false(size(t)); % Boolean array for treatment state

    for i = 1:length(t)
        [~, treatment_status(i), direction(i)] = adp_treat_rate(t(i), r1, r_treat, y(i), v(i), PSA_threshold_low, PSA_threshold_high);
        %fprintf('At time %.2f, PSA = %.2f, Treatment = %d\n', t(i), y(i) + v(i), treatment_status(i));

    end

    figure
    hold on
    plot(t,y,'r', 'DisplayName', 'Sensitive Population')
    plot(t,v,'b', 'DisplayName', 'Resistant Population')
    plot(t,total_count, 'g', 'DisplayName', 'Total')

    % semilogy(t,y,'r', 'DisplayName', 'Sensitive Population')
    % semilogy(t,v,'b', 'DisplayName', 'Resistant Population')
    % semilogy(t,total_count, 'g', 'DisplayName', 'Total')


    %set(gca, 'YScale', 'log')

    x_limits = xlim;
    y_limits = ylim;

    treatment_regions = [0; treatment_status(:); 0];
    edges = diff(treatment_regions);
    start_idx = find(edges == 1);
    end_idx = find(edges == -1); %-1;

    for i = 1:length(start_idx)
        x_start = t(start_idx(i));
        % x_end = t(end_idx(i));

        % Determining x_end to avoid indexing errors:
        if end_idx(i) > length(t)
            x_end = t(end); % Assign the last value of t if index is out of bounds
        else
            x_end = t(end_idx(i)); % Assign the valid index value
        end

        fill([x_start, x_end, x_end, x_start], [0, 0, y_max, y_max], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    end

    treatment_box = fill([0, 1, 1, 0], [-1, -2, -2, -1], [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'on');
    set(treatment_box, 'DisplayName', 'Treatment Applied');



    % if ~isnan(stop_time)
    %     % Progression occurred
    %     plot(stop_time, prog_pnt, 'ro', 'MarkerSize', 8, 'DisplayName', 'Progression Threshold');
    %     legend('Total Tumor Population', 'Progression Threshold', 'Location', 'Best');
    %     title(sprintf('Tumor Progression: Stopped at t = %.1f days', stop_time));
    % else
    %     % No progression
    %     legend('Total Tumor Population', 'Location', 'Best');
    %     title('Tumor Progression: No progression within 600 days');
    % end




    xlabel('Time, (days)')
    ylabel('Cell Count')
    title('Adaptive Treatment Exponential Model');
    legend show;
    ylim([-10,y_max])
    

else
    error('Invalid treatment method specified');
end 