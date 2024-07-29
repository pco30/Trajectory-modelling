function [znext] = stepRungeKutta(t, z, dt, m, start)
% stepRungeKutt computes one step using the RK4 method, which is an
% extension of Eulers method but with 4 derivatives, therefore is more accurate.
% 
% ZNEXT = stepRungeKutta(T,Z,DT,START) computes the state vector ZNEXT at the next
% time step T+DT using stateDeriv

% 4 approximations
A = dt * stateDeriv(z, t, m, start);

B = dt * stateDeriv(z + A/2, t+dt/2, m, start);

C = dt * stateDeriv(z + B/2, t+dt/2, m, start);

D = dt * stateDeriv(z + C, t+dt, m, start);


% Weighted average of approximations to compute next step
znext = z + (A+2*B+2*C+D)/6;

