%exponential model main file
%(I am assuming cell type 1 is sensitive)
r1 = 0.02; % inate growth rate
r2 = 0.015; % inate growth rate
w12 = 0.001;  
w21 = 0.001;  
N1_0 = 10000; % starting population sensitive cells
N2_0 = 10000; % starting pop. resistant cells
PSA_0 = N1_0 + N2_0; % Initial PSA level (sum of initial cell counts)
time_step = 1; % Time step for the simulation (e.g., 1 day)
t = 0:time_step:1000; % Time array

% Parameters for adaptive therapy
PSA_threshold_high = 1.1 * PSA_0; % When treatment is restarted
PSA_threshold_low = 0.9 * PSA_0; % When treatment is withdrawn
delta = 0.02; % Reduction in growth rate due to treatment

% Initial conditions vector:
y0 = [N1_0; N2_0; PSA_0];

% Choose treatment method
treatment_method = 'adaptive'; % options: 'no_treatment', 'continuous', 'metronomic', 'adaptive'

% Solve ODEs based on treatment method
if strcmp(treatment_method, 'no_treatment')
    [T, Y] = ode45(@(t, y) no_treatment(t, y, r1, r2, w12, w21), t, y0);
elseif strcmp(treatment_method, 'continuous')
    [T, Y] = ode45(@(t, y) continuous_treatment(t, y, r1, r2, w12, w21, delta), t, y0);
elseif strcmp(treatment_method, 'metronomic')
    [T, Y] = ode45(@(t, y) metronomic_treatment(t, y, r1, r2, w12, w21, delta), t, y0);
elseif strcmp(treatment_method, 'adaptive')
    [T, Y] = ode45(@(t, y) adaptive_treatment(t, y, r1, r2, w12, w21, PSA_threshold_high, PSA_threshold_low, delta), t, y0);
else
    error('Invalid treatment method specified.');
end

% Plot results
figure;
plot(T, Y(:, 1), 'r', 'DisplayName', 'Cancer Type 1');
hold on;
plot(T, Y(:, 2), 'b', 'DisplayName', 'Cancer Type 2');
plot(T, Y(:, 3), 'g', 'DisplayName', 'PSA Level');
xlabel('Time (days)');
ylabel('Cell Population');
legend;
title(['Tumor Growth with ', treatment_method, ' Treatment']);


%%
% Test setup
y_test = [10; 5; 2]; % Example initial conditions
t_test = 0; % Example time
r1 = 0.1; r2 = 0.2; 
w12 = 0.01; w21 = 0.02;
PSA_threshold_high = 50; PSA_threshold_low = 10; 
delta = 0.05;

% Call function
dydt_test = adaptive_treatment(t_test, y_test, r1, r2, w12, w21, PSA_threshold_high, PSA_threshold_low, delta);

% Display result
disp(dydt_test); % Should display a 3x1 column vector