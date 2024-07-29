function [znext] = stepRungeKutta(t, z, dt)
% stepRungeKutt computes one step using the RK4 method, which is an
% extension of Eulers method but with 4 derivatives, therefore is more accurate.
% 
% ZNEXT = stepRungeKutta(T,Z,DT) computes the state vector ZNEXT at the next
% time step T+DT using stateDeriv

% 4 approximations

A = dt * stateDeriv(z, t); 

B = dt * stateDeriv(z + A/2, t+dt/2);

C = dt * stateDeriv(z + B/2, t+dt/2);

D = dt * stateDeriv(z + C, t+dt);

% Weighted average of approximations to compute next step
znext = z + (A+2*B+2*C+D)/6;

