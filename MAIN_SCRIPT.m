clear all
dbstop if error
format compact

%%% LOAD DATA
load('covid19_case_study')
xdata = day; % temporal data
ydata = moving_av./Nt; % daily case data

%%% SET ENVIRONMENT MODEL PARAMETERS
[trecord] = set_parameters(N,xdata(end));

%%%% INITIAL CONDITIONS
S0 = (N-40-100)/Nt; % susceptible
E0 = 65/Nt; % exposed
I0 = 40/Nt; % infected
R0 = 0; % recovered
C0 = 0; % airborne concentration
D0 = 0; % deposited amount

IC = [C0; D0; S0; E0; I0; R0];

%%% ABC
Naccept = 1;%1500; % number of samples to accept
Nparms = 5; % number of parameters to infer
tol = 1e-8; % tolerance required
discrepany = @(x) 1./length(ydata).*sum((x-ydata).^2); % discrepancy function
accepted = 0; % counter for number of samples accepted
rejected = 0; % counter for number of samples rejected
samples_accepted = zeros(Naccept,Nparms); % store accepted samples
errors = zeros(Naccept,1); % store errors (discrepancy) of accepted samples 

%%% BOUNDS FOR UNIFORM PRIOR DISTRIBUTIONS 
Beta0 = [0.1 0.4]; % units of day^{-1}
Gamma = [1/6 1/3]; % units of day^{-1}
Delta = [1/12 1/10]; % units of day^{-1}
Beta_air = [0.5 5]*1e-7; % units of virus
Beta_sfc = [0.5 5]*1e-7; % units of virus
% NOTE: If using standard SEIR model, set Beta_air = [0 0] and Beta_sfc =
% [0 0]. 

%%% REJECTION SAMPLING
while sum(accepted) < Naccept
    
    %%% EXTRACT PROPOSAL FROM PRIOR
    proposal = get_proposal(Beta0,Beta_air,Beta_sfc,Gamma,Delta);
    
    %%% SOLVE SIMULATION
    if sum(Beta_air + Beta_sfc) == 0
        proposal_days = proposal.*(60*60*24); % solve in days
        [x,x_sol,R0_sol] = solver_SEIR(proposal_days,IC,xdata(2:end),N,Nt);
    else
        [x,x_sol,R0_sol] = solver_SEIRe(proposal,IC,trecord,N,Nt); % solve in seconds
    end
    
    %%% COMPUTE DISCREPANCY
    disc = discrepany(x);
    
    if accepted >= Naccept
        continue
    end
    
    %%% ACCEPT OR REJECT SAMPLE
    if disc <= tol
        accepted = accepted + 1;
        samples_accepted(accepted,:) = proposal;
        errors(accepted) = disc;
    else
        rejected = rejected + 1;
    end
    
    
    if mod(accepted,10) == 0
        disp(accepted)
        save(['N_' num2str(accepted) '_coupled'])
    end
    
    
end

%%% FIND MIN ERROR AND SAVE WORKSPACE
rounded_samples = samples_accepted.*(60*60*24);
[~,idx] = min(errors);
acceptance_rate = accepted./(rejected+accepted);
if sum(Beta_air + Beta_sfc) == 0
    parms_fit = samples_accepted(idx,:).*(60*60*24); % store in days^-1
    [x_fit,x,R0_fit] = solver_SEIR(parms_fit,IC,xdata(2:end),N,Nt);
    error_fit = discrepany(x_fit); % compute discrepancy
    parms_fit(2:3) = 1./parms_fit(2:3); % put into days
    save(['N_' num2str(accepted) '_complete_uncoupled_ma']) % save solution
else
    parms_fit = samples_accepted(idx,:);
    parms_fit(1:3) = samples_accepted(idx,1:3).*(60*60*24); % store in days^-1
    [x_fit,x,R0_fit] = solver_SEIRe(parms_fit./(3600*24),IC,trecord,N,Nt);
    error_fit = discrepany(x_fit); % compute discrepancy
    parms_fit(2:3) = 1./parms_fit(2:3); % put into days
    save(['N_' num2str(accepted) '_complete_coupled_ma']) % save solutions 
end

%%%% PLOT RESULTS
plot_compare_data(x,xdata,ydata.*Nt)
plot_accepted_samples(rounded_samples,parms_fit,errors)

% assign storage for summary statistics 
modes = zeros(Nparms,1);
mus = zeros(Nparms,1);
sigmas = zeros(Nparms,1);
medians = zeros(Nparms,1);

[modes(1),mus(1),sigmas(1),medians(1)] = plot_posterior(rounded_samples(:,1),'\beta_0');
[modes(2),mus(2),sigmas(2),medians(2)] = plot_posterior(1./rounded_samples(:,2),'\gamma^{-1}');
[modes(3),mus(3),sigmas(3),medians(3)] = plot_posterior(1./rounded_samples(:,3),'\delta^{-1}');
[modes(4),mus(4),sigmas(4),medians(4)] = plot_posterior(rounded_samples(:,4),'\beta_{\mathrm{air}}');
[modes(5),mus(5),sigmas(5),medians(5)] = plot_posterior(rounded_samples(:,5),'\beta_{\mathrm{sfc}}');

%%% PRINT RESULTS TO SCREEN
fprintf('\nFitted parameters: %f %f %f %f %f\n', parms_fit')
fprintf('Fitted error: %f\n', error_fit)
fprintf('Acceptance rate: %f\n', acceptance_rate)
fprintf('Posterior modes: %f %f %f %f %f\n', modes')


