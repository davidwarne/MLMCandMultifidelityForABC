function [FP, FN, TP, TN] = MultifidelityROC(N,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx)

%% Early accept/reject multifidelity for approximate Bayesian computaion to compute
% the posterior mean
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
%
% Outputs:
%
% Authors:
%   Thomas P. Prescott[1] (prescott@maths.ox.ac.uk)
%   David J. Warne[2,3,4] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] Mathematical Institute, University of Oxford, UK
%   [2] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [3] Centre for Data Science, Queensland University of Technology, Autralia
%   [4] ARC Centre of Excellence for Mathematical and Statistical Frontiers

% initialise
% theta = [];
FP = 0;
FN = 0;
TP = 0;
TN = 0;

for i = 1:N
    % generate trial from the prior
    theta_trial = p();
    % simulate approximate data using these parameters
    D_s_approx = s_approx(theta_trial);
    dist_approx = rho_approx(D_s_approx);

    % simulate exact model
    D_s = s(theta_trial);
    dist = rho(D_s);
%    theta = [theta,theta_trial];

    if dist_approx<epsilon_approx
        if dist<epsilon
            TP = TP+1;
        else
            FP = FP+1;
        end
    else
        if dist<epsilon
            FN = FN+1;
        else
            TN = TN+1;
        end
    end
end

FP = FP/N;
FN = FN/N;
TP = TP/N;
TN = TN/N;

end