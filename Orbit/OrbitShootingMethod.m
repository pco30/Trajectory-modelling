function  [start,z,t] = OrbitShootingMethod(alt)
%   ShootingMethod is a Boundary Value Problem (BVP) solver which 
%   unsuccesfullly attempts to solve the altitude at which thrust begins


start(1) = 10000; % First altitude guess (m)
start(2) = 186050; % Seconed altitude guess (m)

%% Manually obtaining third guess for start
[alt1] = ivpSolver(0,54.0328,2500,0.1,start(1)); % Solving IVP of first guess
z1 = max(alt1(5,:)); % Max altitude for first guess

[alt2] = ivpSolver(0,54.0328,2500,0.1,start(2)); % Solving IVP of second guess
z2 = max(alt2(5,:)); % Max altitude for second guess

e1 = z1-alt; % Error of first guess
e2 = z2-alt; % Error of second guess

%Calculating equation of the error-angle line to calculate new guess for
%start
m = (e2-e1) / (start(2)-start(1)); % Gradient
c = e2 - m*start(2); % Y-intercept
start(3) = -c / m; % New start guess

%% Repeating the steps above, using the the nth and n-1th start term to
% calculate the next start term. Loop terminates when start is 99.99%
% accurate.
n=2;
while abs(e2) > alt*0.0001

    [alt1] = ivpSolver(0,54.0328,2500,0.1,start(n)); 
    z1 = max(alt1(5,:));
    
    [alt2] = ivpSolver(0,54.0328,2500,0.1,start(n+1));
    z2 = max(alt2(5,:));
    
    e1 = z1-alt;
    e2 = z2-alt;
    
    m = (e2-e1) / (start(n+1)-start(n));
    c = e2 - m*start(n+1);
    
    start(n+2) = -c / m; % New start guess

    n = n + 1;

end

[z,t]=ivpSolver(0,54.0328,v,0.1,start(end));


%% Plotting results [as code loops infintely as it results in NaN, this stage never occurs]
xlim([z(1,1) z(1,end-1)+10^4])
ylim([z(3,end-1)-10^4 max(z(3,:))+10^4])
hold on
r = 6.3878*10^6;
plot(r*cos(0:2*pi/1000:2*pi),r*sin(0:2*pi/1000:2*pi)-r)