function [E,V,ESS] = ABCMultifidelity(N,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx,eta1,eta2,f)
%% Early accept/reject multifidelity for approximate Bayesian computation 
% to compute expectation E[f(theta)] with respect to the ABC posterior
%
% Inputs:
%    N - Number of samples
%    p - prior distribution sampler, 
%    s - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - discrepancy acceptance threshold
%    s_approx - function that generates simulated approximate data give a parameters set
%    rho_approx - approximate discrepancy metric, treated as a function of simulated data only
%    epsilon_approx - approximate discrepancy treshold
%    eta1 - continuation probability given acceped approximate sample
%    et2 - continuation probability given rejected approximate sample
%    f - functional to compute E[f(theta)] w.r.t to ABC posterior measure
%
% Outputs:
%    E - Monte Carlo estimate of E[f(theta)]
%    V - Estimator variance (based on delta method approximation)
%    ESS - Effective Sample Size
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

% initialise
theta = [];
w = [];
for i = 1:N
    % generate trial from the prior
    theta_trial = p();
    % simulate approximate data using these parameters
    [D_s_approx, couple_arg_1, couple_arg_2, couple_arg_3] = s_approx(theta_trial);
    % compute early accept/reject weight
    w(i) = (rho_approx(D_s_approx) <= epsilon_approx);
    % compute continuation probility
    eta =  eta1*w(i) + eta2*(1- w(i));
    % continue with probability eta
    if unifrnd(0,1) < eta
        % simulate exact model
        D_s = s(theta_trial, couple_arg_1, couple_arg_2, couple_arg_3);
        % update weights
        w(i) = w(i) + ((rho(D_s) <= epsilon) - w(i))/eta;
    end
    theta = [theta,theta_trial];
end
F = f(theta);
% compute Multifiedlity estimator
E = sum(w.*F)/(sum(w));
% compute Variance applroximation
mu_w = mean(w);
mu_wf = mean(w.*F);
c = cov(w,w.*F);
V = (var(w)*(mu_wf/mu_w)^2 + var(w.*F) - 2*c(1,2)*(mu_wf/mu_w))/(N*mu_w^2);
% effective sample size (proportional to ESS \propto 1/V)
ESS = (sum(w)^2)/sum(w.*w);

