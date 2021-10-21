function [theta] = ABCRejectionSampler(N,p,s,rho,epsilon)
%% Rejection Sampler for approximate Bayesian computaion
%
% Inputs:
%    N - the number of ABC posterior samples
%    p - function that generates iid samples from the parameter joint
%        prior distribution
%    s - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - the discrepancy acceptance threshold
%
% Outputs:
%    theta - a matrix of ABC posterior samples
%
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers

theta = [];
while size(theta,2) < N
    % generate trial from the prior
    theta_trial = p();
    % generate simulated data using these parameters
    D_s = s(theta_trial);
    % accept/reject
    if rho(D_s) <= epsilon
        theta = [theta,theta_trial];
    end
end
