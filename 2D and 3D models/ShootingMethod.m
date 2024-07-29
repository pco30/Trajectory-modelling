function  [theta,z,t] = ShootingMethod(alt,v,dt)
% ShootingMethod is a numerical root finder for the boundary value problem 
% (BVP) of the desired apogee of the projectile, alt.

theta(1) = 0.1; % First angle guess, in degrees (°)
theta(2) = 60; % Seconed angle guess, in degrees (°)

%% Manually obtaining third guess for theta
[alt1] = ivpSolver(0,theta(1),v,dt); % Solving IVP of first guess
z1 = max(alt1(5,:));% Max value of altitude of first theta guess

[alt2] = ivpSolver(0,theta(2),v,dt); % Solving IVP of second guess
z2 = max(alt2(5,:)); % Max value of altitude of second theta guess

e1 = z1-alt; % First guess error
e2 = z2-alt; % Second guess error

% Calculating equation of the error-angle line to calculate new guess for
% theta
m = (e2-e1) / (theta(2)-theta(1)); % Gradient
c = e2 - m*theta(2); % Y-intercept
theta(3) = -c / m; % New theta guess

%% Repeating the steps above, using the the nth and n-1th theta term to
% calculate the next theta term. Loop terminates when theta is 99.99999%
% accurate.
n=2;
while abs(e2) > alt*0.0000001

    [alt1] = ivpSolver(0,theta(n),v,dt); 
    z1 = max(alt1(5,:));
    
    [alt2] = ivpSolver(0,theta(n+1),v,dt);
    z2 = max(alt2(5,:));
    
    e1 = z1-alt;
    e2 = z2-alt;
    
    m = (e2-e1) / (theta(n+1)-theta(n));
    c = e2 - m*theta(n+1);
    
    theta(n+2) = -c / m; % New theta guess

    n = n + 1;

end

[z,t]=ivpSolver(0,theta(end),v,dt); % Finding state vector, z, for correct theta