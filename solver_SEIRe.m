function [I_all,x,R0] = solver_SEIRe(parms,IC,trecord,N,Nt)
% FUNCTION SOLVER_SEIRe
%
% solves the SEIR-e model using Euler's Method
% parms     === parameters
% IC        === initial condition 
% trecord   === time to record solution at (record daily)
% N         === number of people 
% Nt        === total population of Victoria

global dt t_end
global p
global phi psi kappa_air kappa_sfc epsilon nu
global N_air velocity_dep h

%%% EXTRACT PARAMETERS
Beta = parms(1);
Gamma = parms(2);
Delta = parms(3);
Beta1 = parms(4);
Beta2 = parms(5);

load('temp_data') % load temperature data
tsolve = 0:dt:t_end; % time array to solve over
dtau = (tsolve(2)-tsolve(1))*Delta; % non-dimensional time step 
counter = 1; % counter for solution storage

%%% COMPUTE PARAMETERS IN ENVIRONMENT MODEL
eta = lambda_surface_fit(1);
gamma = lambda_air_fit(1) + N_air + velocity_dep./h;
Phi = gamma + phi;
Psi = eta + psi;

%%% EXTRACT INITIAL CONDITIONS 
C_old = IC(1); D_old = IC(2); 
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
R0(1) = (Beta)/(Delta*N) + Beta1.*kappa_air/Delta*sum(nu./Phi.*p) + Beta2.*kappa_sfc./(Psi.*Delta).*sum(epsilon.*nu./Phi);

%%% TIME LOOP
for tidx = 2:length(tsolve)
    t = tsolve(tidx); % real time in seconds 
    
    idx = ceil(t/(3600*24)); % rounded real time in days
    
    %%% UPDATE PARAMETERS IN ENVIRONMENT MODEL
    Phi = lambda_air_fit(idx) + N_air + velocity_dep./h + phi;
    Psi = lambda_surface_fit(idx) + psi;
    
    %%% SOLVE ODE USING EULER METHOD
    C = C_old + dtau*(-Phi/Delta.*C_old + I_old); C_path = sum(C.*nu.*p); % number of pathogens
    D = D_old + dtau*(C_old - Psi/Delta.*D_old); D_path = sum(D.*epsilon.*nu); % number of pathogens    
    Susceptible = S_old + dtau*(-Beta/Delta*S_old*I_old*Nt/N ...
        - Beta1/Delta^2*kappa_air*Nt*S_old*C_path ...
        - Beta2/Delta^3*kappa_sfc*Nt*S_old*D_path);
    Exposed = E_old + dtau*(Beta/Delta*S_old*I_old*Nt/N ...
        + Beta1/Delta^2*kappa_air*Nt*S_old*C_path ...
        +  Beta2/Delta^3*kappa_sfc*Nt*S_old*D_path ...
        - Gamma/Delta*E_old);
    Infected = I_old + dtau*(Gamma/Delta*E_old - I_old) ;
    Recovered = R_old + dtau*(I_old);
    
     %%% UPDATE SOLUTION
    S_old = Susceptible;
    I_old = Infected;
    R_old = Recovered;
    E_old = Exposed;
    C_old = C;
    D_old = D;
    
    %%% STORE SOLUTIONS
    if abs(trecord(counter)-t) <= 1e-3
        S_all(counter+1) = Susceptible;
        I_all(counter+1) = Infected;
        R_all(counter+1) = Recovered;
        E_all(counter+1) = Exposed;
        
        Chat = nu.*Nt./Delta;
        Dhat = epsilon.*nu.*Nt./Delta^2;
        C_all(counter+1) = sum(C.*Chat);
        D_all(counter+1) = sum(D.*Dhat);
        
        R0(counter+1) = (Beta)/(Delta*N) + Beta1.*kappa_air/Delta*sum(nu./Phi.*p) ...
        + Beta2.*kappa_sfc./(Psi.*Delta).*sum(epsilon.*nu./Phi); % compute R0
    
        fprintf(['\n beta_1 = ' num2str(Beta1) '; t = ' num2str(trecord(counter)/(24*3600)) ' days\n'])
        
        counter = counter + 1;
    end
    
end

x = [S_all.*Nt E_all.*Nt I_all.*Nt R_all.*Nt C_all D_all]; % combine solutions to extract



end