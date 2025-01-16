%Trying to investigate how the equilibria of the system of dif. eqns.
%relies on the two phenotypic switching rates. Setting r1 and r2 to be
%constant for this investigation

%initialising intrinsic rates to 1 for the purpose of investigating
%equilibria
r1_eq = 1;
r2_eq = 1;


%initialising omega values
w_start = 0;
w_end = 1;
w_int = 100;

w12_vals = linspace(w_start, w_end, w_int);
w21_vals = linspace(w_start, w_end, w_int);

N1_eq_vals = zeros(length(w12_vals), length(w21_vals));
N2_eq_vals = zeros(length(w12_vals), length(w21_vals));
stability_eval = zeros(length(w12_vals), length(w21_vals));


for i = 1:length(w12_vals)
    for j = 1:length(w21_vals)
        %Assigning w12 and w21 for each iteration
        w12 = w12_vals(i);
        w21 = w21_vals(j);

        %creating a vector for the system of differential equations dN1/dt
        %and dN2/dt to solve for possible equilibria (when both dN1/dt and
        %dN1/dt = 0)

        vec_system = @(N) [
            r1_eq*N(1) - w12*N(1) + w21*N(2);
            r2_eq*N(2) - w21*N(2) + w12*N(1)];

        [Eq_vals, ~, exitflag] = fsolve(vec_system, [100, 100]); 
        %using fsolve to find roots of system - initial N1 and N2 guess for
        %iterations [a, b]

        if exitflag > 0 %exitflag tests the reason fsolve stopped. if >0 means a solution was found
            N1_eq_vals(i,j) = Eq_vals(1);
            N2_eq_vals(i,j) = Eq_vals(2);

            %Stability analysis of equilibria if equilibria is found:
            
            %Creating Jacobean matrix
            J = [r1_eq - w12, w21;
                w12, r2_eq - w21]; 
            %Determining eigen-values of Jacobean matrix
            eigenvals = eig(J);

            if all(real(eigenvals) < 0) %stable
                stability_eval(i,j) = 1;
            else
                stability_eval(i,j) = -1; %unstable
            end


        else
            N1_eq_vals(i,j) = NaN;
            N2_eq_vals(i,j) = NaN;
        end

    end
end


% Visualize results
figure;
contourf(w12_vals, w21_vals, N1_eq_vals'); % Equilibria for N1
xlabel('\omega_{12}');
ylabel('\omega_{21}');
title('Equilibria for N1');
colorbar;


figure;
contourf(w12_vals, w21_vals, stability_eval);
xlabel('\omega_{12}');
ylabel('\omega_{21}');
title('Equilibria for N1');
colorbar;
