function [dz] = stateDeriv(z, t, m, start)
% Calculate the state derivatives for the projectile

% 
%     DZ = stateDeriv(T,Z) computes the derivative DZ = [Vx; Ax; Vy; Ay] of the 
%     state vector Z = [Dx; Vx; Dy; Vy], where D displacement, V is velocity,
%     and A is acceleration for x and y.

% Variable definitions
G = 6.674*10^-11; % Gravitational Constant (m^3/kg/s^2)
M = 5.972*10^24; % Mass of Earth (Kg)
r = 6.3878*10^6; % Radius of Earth (m)


% Calculations
alt = sqrt(z(1,end)^2+(z(3,end)+r)^2)-r; % Altitude (m)

[rho,T,p] = atmosEarth(alt); % Air density as a function of altitude (kg/m^3)

F = thrust(200,3000,150000,p,0.75^2*pi); % Thrust force produced

thetaA=atan2(z(3,end)+r,z(1,end)); % Angle relative to Earth's centre (rad)

thetaR=atan2(z(4,end),z(2,end)); % Orientation angle of projectile (rad)

v = sqrt(z(4,end)^2+z(2,end)^2); % Apsolute velocity (m/s)

grav = -(G*M) / ((r+alt)^2); % Gravity as a function of altitude (m/s^2)


% Asigning velocity to the state derivative
dzVx = z(2); 
dzVy = z(4);

% Calculating acceleration

if  and(alt>start,m>4500) % Projectile motion ODEs with thrust
    dzAy = grav * sin(thetaA)  - 1 / (2*m) * rho * 0.1 * pi * v^2*sin(thetaR) + F/m * sin(thetaR);
    dzAx = grav * cos(thetaA)  - 1 / (2*m) * rho * 0.1 * pi * v^2*cos(thetaR) + F/m * cos(thetaR);

else % Unpowered projectile ODEs
    dzAy = grav * sin(thetaA) - 1 / (2*m) * rho * 0.1 * pi * v^2*sin(thetaR);
    dzAx = grav * cos(thetaA) - 1 / (2*m) * rho * 0.1 * pi * v^2*cos(thetaR);
end

% State derivative vector%
dz = [dzVx; dzAx; dzVy; dzAy];