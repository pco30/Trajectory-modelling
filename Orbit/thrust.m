function F = thrust(dm,ve,pe,p0,A)
%% thrust  Model for rocket thrust
%
% [F] = thrust(DM,VE,PE,P0,A) outputs the modelled values of rocket thrust 
% F in Newtons, given the input parameters ME for mass flow rate in kg/s, 
% VE for exhaust velocity in m/s, PE for exit pressure in Pa, P0 for 
% external pressure in Pa, and A for exhaust cross-sectional area in m^3.
%
% Model from "Spacecraft Systems Design and Engineering", Encyclopedia of 
% Physical Science and Technology
% 
% Alan Hunter & Lonox Jianqiang Huang, University of Bath, Oct 2022.


%% Momentum thrust
Fm = dm * ve;


%% Pressure thrust
Fp = (pe-p0) * A;


%% Total thrust
F = Fm + Fp;