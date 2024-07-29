function [z,t] = ivpSolver(t0,theta,v,dt,start)
% ivpSolver    Solve an initial value problem (IVP) and plot the result
% 
%     [T,Z] = ivpSolver(T0,THETA,V,DT,START) computes the IVP solution using a step 
%     size DT, beginning at time T0 and initial velocity at angle THETA and
%     thrust beginning at altitude START
%     The solution is output as a time vector T and a matrix of state 
%     vectors Z.


% Z ROWS
% Row 1 - X displacement (m)
% Row 2 - X velocity (m/s)
% Row 3 - Y displacement (m)
% Row 4 - Y velocity (m/s)
% Row 5 - Altitude (m)
% Row 6 - Mass of rocket (Kg)

% Set initial condition
t(1) = t0;
z(:,1) = [0;v*cosd(theta);0;v*sind(theta);0];
z(6,1) = 25000;


r = 6.3878*10^6; % Earth Radius

% Continue stepping until projectile has hit the ground (0m altiude), or
% 8000s has passed
n=1;
while and(t(end)<8000,z(5,end)>=0)
    % Increment the time vector by one time step
    t(n+1) = t(n) + dt;
    
    % Apply Runge Kutta method for one time step
    [z(1:4,n+1)] = stepRungeKutta(t(n), z(1:4,n), dt, z(6,end),start);
    
    z(5,n+1) = sqrt(z(1,end)^2+(z(3,end)+r)^2)-r; % Calculation of Altitude
    
    % Projectile Mass calcultations

    % Locating time at which thrust begins
    a = z(5,:)>start;
    b = diff(a);
    [c,d] = max(b);

    %
    if and(z(5,end)>start,z(6,n)>4500)
        z(6,n+1) = 25000 - 200*(t(n)-t(d)); % Thrust stage burning mass

    elseif z(5,end)<start
        z(6,n+1) = 25000;  % Constant mass up to START

    else
        z(6,n+1) = z(6,n); % Constant mass after thrust of 4500Kg

    end

    n = n+1; % Setting next loop
end


%% Plotting results
plot(z(1,:),z(3,:),'r')
hold on
plot(r*cos(0:2*pi/1000:2*pi),r*sin(0:2*pi/1000:2*pi)-r,'g')
a=r+2000000;
plot(a*cos(0:2*pi/1000:2*pi),a*sin(0:2*pi/1000:2*pi)-r,'k',LineStyle='--')
axis equal
legend('Projectile Path','Earth','2000Km Orbit')
xlabel('X-Displacement (m)')
ylabel('Y-Displacement (m)')