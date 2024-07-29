function [z,t] = ivpSolver(t0,theta,v,dt)
% ivpSolver    Solve an initial value problem (IVP) and plot the result
% 
%     [T,Z] = ivpSolver(T0,THETA,V,DT) computes the IVP solution using a step 
%     size DT, beginning at time T0 and initial velocity at angle THETA.
%     The solution is output as a time vector T and a matrix of state 
%     vectors Z.

% Z ROWS
% Row 1 - X displacement (m)
% Row 2 - X velocity (m/s)
% Row 3 - Y displacement (m)
% Row 4 - Y velocity (m/s)
% Row 5 - Altitude (m)

% Set initial conditions
t(1) = t0;
z(:,1) = [0;v*cosd(theta);0;v*sind(theta);0];


r = 6.3878*10^6; % Earth Radius (m)

n=1;

% Continue stepping until projectile has hit the ground (0m altiude)
while z(5,end)>=0
    % Increment the time vector by one time step
    t(n+1) = t(n) + dt;
    
    % Apply Runge Kutta method for one time step
    [z(1:4,n+1)] = stepRungeKutta(t(n), z(1:4,n), dt);
    
    z(5,n+1) = sqrt(z(1,end)^2+(z(3,end)+r)^2)-r; % Done here because rungekutta was creating error for altitude
    
    n = n+1;
end
