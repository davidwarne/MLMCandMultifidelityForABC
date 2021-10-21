function [Y_obs] = GenerateCoupledObservations(bcrn,k,X0,N,I,t,sigma,D_obs,P_obs,zeta_obs)
%% Generate Observations
% Generates N Observations of state dimensions I of stochastic process X(t) conditional
% on approximate observations with rate parameters k and initial condition X(0) = X0 = X0. 
% Observations are taken at time points t(1),..,t(Nt) with observation error 
% Y_obs(t(i)) = X(t(i)) + zeta, where zeta ~ N(0,sigma^2)
%
% Inputs:
%    bcrn - a BCRN structure
%    k    - true kinetic rate parameters for the model
%    X0   - true initial condition of the BCRN
%    N    - number of sample paths in dataset
%    I    - vector of observable chemical species indices
%    t    - vector of observation times
%    sigma - standard deviation of observation error
%    D_obs - the coarse-grained (tau-leap) reaction-time interval lengths
%    P_obs - the Poisson number of reaction firings in each interval
%    zeta_obs - the Gaussian measurement noise
%
%    Latter three inputs are the three optional outputs of:
%    [Y_obs, D_obs, P_obs, zeta_obs] = GenerateApproxObservations(bcrn,k,X0,N,I,t,sigma,tau)
%
% 
% Outputs:
%     a Table of observations Y_obs
%
% Authors:
%   Thomas P. Prescott[1] (tprescott@turing.ac.uk)
%   David J. Warne[2,3,4] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] The Alan Turing Institute, London, UK
%   [2] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [3] Centre for Data Science, Queensland University of Technology, Autralia
%   [4] ARC Centre of Excellence for Mathematical and Statistical Frontiers

Nt = length(t);
Y_obs = zeros(length(I),Nt,N);
bcrn.k = k;
for i=1:N
    % stochastic simulation
    [X_r,t_r] = CoupledNextReactionMethod(bcrn,t(Nt),D_obs(:,:,i),P_obs(:,:,i));
    % just store the state of selected species and times
    for j=1:Nt
        [J] = find(t_r <= t(j));
        Y_obs(:,j,i) = X_r(I,J(end));
    end
end
% perturb with noise
% zeta = normrnd(0,sigma,size(Y_obs));
Y_obs = Y_obs + zeta_obs(:,:,1:N);
