% Parameters (example values)
r1 = 1; r2 = 0.5;
omega12 = 0.3; omega21 = 0.2;

% Time span and initial conditions
tspan = [0 50];
N0 = [1; 1]; % Initial populations

% Define the system of ODEs
dNdt = @(t, N) [r1*N(1) - omega12*N(1) + omega21*N(2);
                 r2*N(2) - omega21*N(2) + omega12*N(1)];

% Solve the system using ode45
[t, N] = ode45(dNdt, tspan, N0);

% Total population
Y = sum(N, 2);

% Plot results
figure;
plot(t, N(:,1), '-r', 'DisplayName', 'N1');
hold on;
plot(t, N(:,2), '-b', 'DisplayName', 'N2');
plot(t, Y, '-k', 'DisplayName', 'Total Population');
xlabel('Time');
ylabel('Population');
legend;
title('Population Dynamics');
