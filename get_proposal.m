function out = get_proposal(Beta0,Beta_air,Beta_sfc,Gamma,Delta)
% FUNCTION GET_PROPOSAL
%
% extract randomly generate samples from uniform distribution 
% Beta0     === bounds of uniform distribution for Beta0
% Beta_air  === bounds of uniform distribution for Beta_air
% Beta_sfc  === bounds of uniform distribution for Beta_sfc
% Gamma     === bounds of uniform distribution for Gamma
% Delta     === bounds of uniform distribution for Delta

beta = Beta0(1) + (Beta0(2)-Beta0(1))*rand;
gamma = Gamma(1) + (Gamma(2)-Gamma(1))*rand;
delta = Delta(1) + (Delta(2)-Delta(1))*rand;
beta_air = Beta_air(1) + (Beta_air(2)-Beta_air(1))*rand;
beta_sfc = Beta_sfc(1) + (Beta_sfc(2)-Beta_sfc(1))*rand;

out = [beta/(60*60*24); gamma/(60*60*24); delta/(60*60*24);...
       beta_air; beta_sfc]; 





end