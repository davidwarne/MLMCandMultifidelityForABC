function [eta1,eta2,p_fp,p_fn,p_tp,p_tn,ctilde,cp,cn,rho_dist,rho_dist_approx] = MultifidelityROC(N,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx)
%% Burn-in run to compute false/true positive/negative rates
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
%    eta1 - optimal continuation probability for rho_approx < epsilon_approx
%    eta2 - optimal continuation probability for rho_approx > epsilon_approx
%    p_fp - proportion with (rho_approx < epsilon_approx) & (rho > epsilon)
%    p_fn - proportion with (rho_approx > epsilon_approx) & (rho < epsilon)
%    p_tp - proportion with (rho_approx < epsilon_approx) & (rho < epsilon)
%    p_tn - proportion with (rho_approx > epsilon_approx) & (rho > epsilon)
%    ctilde - expected approximate simulation cost
%    cp - expected exact simulation cost for continuation of positive test
%    cn - expected exact simulation cost for continuation of negative test
%
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
FP = 0;
FN = 0;
TP = 0;
TN = 0;
cp = 0;
cn = 0;
ctilde = 0;
rho_dist = [];
rho_dist_approx = [];
for i = 1:N
    % generate trial from the prior
    theta_trial = p();
    % simulate approximate data using these parameters and time
    tic;
    D_s_approx = s_approx(theta_trial);
    ctilde = ctilde + toc;
    dist_approx = rho_approx(D_s_approx);

    rho_dist_approx = [rho_dist_approx,dist_approx];
    % simulate exact model and time
    tic;
    D_s = s(theta_trial);
    ctemp = toc;
    dist = rho(D_s);

    rho_dist = [rho_dist,dist];
     % compute ROC and expected compute times
    if dist_approx<epsilon_approx
        cp = cp + ctemp;
        if dist<epsilon
            TP = TP+1;
        else
            FP = FP+1;
        end
    else
        cn = cn + ctemp;
        if dist<epsilon
            FN = FN+1;
        else
            TN = TN+1;
        end
    end
end

% ROC estimates 
p_fp = FP/N;
p_fn = FN/N;
p_tp = TP/N;
p_tn = TN/N;

% compuational cost estimates
ctilde = ctilde/N;
cp = (cp/(TP+FP))*(p_tp + p_fp);
cn = (cn/(TN+FN))*(p_tn + p_fn);

% opitmal continuation probabilities (Lemma 4.2 in Prescott and Baker (2020))
R_p = p_fp/(cp/ctilde);
R_n = p_fn/(cn/ctilde);
R_0 = p_tp - p_fp;
eta1 = sqrt(R_p/R_0);
eta2 = sqrt(R_n/R_0);
end
