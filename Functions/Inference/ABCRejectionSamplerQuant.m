function [theta,r] = ABCRejectionSampler(N,p,s,rho,q)
%% Rejection Sampler for approximate Bayesian computaion
%
% Inputs:
%    N - the number of ABC posterior sample trials
%    p - function that generates iid samples from the parameter joint
%        prior distribution
%    s - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    q - quantile for final output
%
% Outputs:
%    theta - a matrix of ABC posterior samples
%    r     - a sorted list of discrepancies r(end) is the effective epsilon
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

theta = [];
r = [];
while size(theta,2) < N
    % generate trial from the prior
    theta_trial = p();
    % generate simulated data using these parameters
    D_s = s(theta_trial);
    theta = [theta,theta_trial];
    r = [r,rho(D_s)];
end
% accept/reject based on quantile
[r, I] = sort(r,'ascend');
theta = theta(:,I);
epsilon = quantile(r,q);
J = find(r > epsilon);
r(J) = [];
theta(:,J) = [];

