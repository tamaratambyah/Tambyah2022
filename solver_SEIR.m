function [I_all,x,R0] = solver_SEIR(parms,IC,trecord,N,Nt)
% FUNCTION SOLVER_SEIR
%
% solves the SEIR model using Euler's Method
% parms     === parameters
% IC        === initial condition 
% trecord   === time to record solution at (record daily)
% N         === number of people 
% Nt        === total population of Victoria


%%% EXTRACT PARAMETERS
Beta = parms(1);
Gamma = parms(2);
Delta = parms(3);

load('temp_data') % load temperature data
tsolve = 0:0.01:trecord(end); % time array to solve over
dtau = (tsolve(2)-tsolve(1))*Delta; % non-dimensional time step
r = Nt/N; % non-dimensional scaling 
counter = 1; % counter for solution storage

%%% EXTRACT INITIAL CONDITIONS 
S_old = IC(3); E_old = IC(4); I_old = IC(5); R_old = IC(6);

%%% STORE SOLUTIONS DAILY
S_all = zeros(length(trecord)+1,1);
E_all = zeros(length(trecord)+1,1);
I_all = zeros(length(trecord)+1,1);
R_all = zeros(length(trecord)+1,1);
C_all = zeros(length(trecord)+1,1); 
D_all = zeros(length(trecord)+1,1);
R0 = zeros(length(trecord)+1,1);

%%% STORE IC 
S_all(:,1) = S_old; 
I_all(:,1) = I_old;
R_all(:,1) = R_old;
E_all(:,1) = E_old;
R0(1) = (Beta)/Delta;

%%% TIME LOOP
for tidx = 2:length(tsolve)
    t = tsolve(tidx); % extract time
    
    %%% SOLVE ODE USING EULER METHOD
    Susceptible = S_old + dtau*(-Beta/Delta*S_old*I_old*r );
    Exposed = E_old + dtau*(Beta/Delta*S_old*I_old*r  - Gamma/Delta*E_old);
    Infected = I_old + dtau*(Gamma/Delta*E_old - I_old) ;
    Recovered = R_old + dtau*(I_old);
    
    %%% UPDATE SOLUTION
    S_old = Susceptible;
    I_old = Infected;
    R_old = Recovered;
    E_old = Exposed;
    
    
    %%% STORE SOLUTIONS
    if abs(trecord(counter)-t) <= 1e-3
        S_all(counter+1) = Susceptible;
        I_all(counter+1) = Infected;
        R_all(counter+1) = Recovered;
        E_all(counter+1) = Exposed;
        R0(counter+1) = (Beta)/Delta;
        counter = counter + 1;
    end
    
end

x = [S_all.*Nt E_all.*Nt I_all.*Nt R_all.*Nt C_all D_all]; % combine solutions to extract


end