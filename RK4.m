function [t,y,v,h] = RK4(f1, f2, a, b, n, N1_0, N2_0)
%function [t,y,v,h,stop_time] = RK4(f1, f2, a, b, n, N1_0, N2_0, prog_pnt)
% RK4 performs the Range-Kutte 4 Method to estimate y and v from two ODEs of the
% form dy/dt = f1(t,y,v) and dv/dt = f2(t,y,v) 
% INPUTS:
%   f1 is RHS of the first ODE 
%   f2 is RHS of the second ODE
%   a and b provide the range for the time interval t, i.e. a<t<b
%   n is the number of steps to get from the a to b 
%   alpha is the initial value of y 
%   beta is the initial value of v 
% OUTPUTS: 
%   t is an array containing the time at each point
%   y is an array containing population of sensitive population at each
%   point
%   v is an array containing population of resistant population at each
%   point

% Form the t array
h = (b-a)/n;    % calculate h, the interval between points
t = a:h:b;

% Initialise the y and v arrays
y = zeros(size(t));
y(1) = N1_0;
v = zeros(size(t));
v(1) = N2_0;

% Perform iterations to fill the y and v arrays
for i = 1:n

    k1 = h*f1(t(i), y(i), v(i));
    m1 = h*f2(t(i), y(i), v(i));

    k2 = h*f1(t(i) + h/2, y(i) + k1/2, v(i) + m1/2);
    m2 = h*f2(t(i) + h/2, y(i) + k1/2, v(i) + m1/2);

    k3 = h*f1(t(i) + h/2, y(i) + k2/2, v(i) + m2/2);
    m3 = h*f2(t(i) + h/2, y(i) + k2/2, v(i) + m2/2);

    k4 = h*f1(t(i) + h, y(i) + k3, v(i) + m3);
    m4 = h*f2(t(i) + h, y(i) + k3, v(i) + m3);

    y(i+1) = y(i) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    v(i+1) = v(i) + 1/6 * (m1 + 2*m2 + 2*m3 + m4);    


    % if y(i+1) + v(i+1) > prog_pnt
    %     stop_time = t(i+1);
    %     t = t(1:i+1);
    %     y = y(1:i+1);
    %     v = v(1:i+1);
    %     break;
    % end
end

end
