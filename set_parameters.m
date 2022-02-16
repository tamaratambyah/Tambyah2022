function [trecord] = set_parameters(N,Ndays)
% FUNCTION SET_PARAMETERS
% 
% assigns the default parameter values stated in DST-Group-RR-XXXX. 
% NOTE: The production rate (e_i) values are assigned from a pre-saved
% workspace. Use the accompany excel spreadsheet to vary the data if
% needed.
% NOTE: Time units is seconds, length units is meters, [] is used to
% denonte units of each parameter

global dt t_end
global p
global phi psi kappa_air kappa_sfc epsilon nu
global N_air velocity_dep h

load('emission_data.mat') % load pre-saved emission data 

%%% TIME PARAMETERS
dt = 5; % time step [seconds]
t_end = Ndays*(24*3600); % final time [seconds] 
trecord = [1:Ndays].*24*3600; % times of record solution at [seconds]

%%% ROOM DIMENSION PARAMETERS
V = 240; % volume of indoor environment [meters^3]
h = 3; % height of indoor environment [meters]

%%% DEPOSITION PARAMETERS
v_max = 0.6; % maximum deposition velocity for conservation [meters/second]
velocity_dep(velocity_dep>v_max) = v_max; % enforce maximum 

%%% SURFACE PICKUP PARAMETERS
A = 10; % area of horizontal surfaces [meters^2]
A_sfc = 0.02; % area of hands (pickup) [meters^2]
R_sfc = 1./(3600); % rate of surface contact [1/seconds]
q = 0.1; % pickup efficiency []
r = 0.34; % transfer efficieny []

%%% INHALATION PARAMETERS
BR = 10*1e-3/(60); % breathing rate [meters^3/second];
threshold = 30; % diameter of particles small enough to remain aloft [microns]
N_air = 2/(3600); % air exchange rate [1/seconds]

%%% OTHER PARAMETERS
alpha = 0.5; % size correction factor @ RH = 40%
effective_diameter = alpha.*diameter; % effective droplet diameter
p = effective_diameter < threshold;  % determine particles less than threshold

%%% COMPUTE ENVIRONMENT MODEL PARAMETERS (time-independent parameters)
nu = emission./V;
epsilon = velocity_dep;
kappa_air = BR;
kappa_sfc = r*q*A_sfc*R_sfc;
phi = p.*BR.*N./V;
psi = q.*A_sfc.*R_sfc.*N./A;

end