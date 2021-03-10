function [E,V,ESS,eta1,eta2] = ABCAdaptiveMultifidelity(N,M,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx,f)
%% Adaptive early accept/reject multifidelity for approximate Bayesian computation
% The function adaptively updates the optimal continuation probabililties based on
% increasingly accurate estimates of ROC properties.
%
% Inputs:
%    N - Number of samples
%    M - Number of burnin samples
%    p - prior distribution sampler, 
%    s - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - discrepancy acceptance threshold
%    s_approx - function that generates simulated approximate data give a parameters set
%    rho_approx - approximate discrepancy metric, treated as a function of simulated data only
%    epsilon_approx - approximate discrepancy treshold
%    f - functional to compute E[f(theta)] w.r.t to ABC posterior measure
%
% Outputs:
%    E - Monte Carlo estimate of E[f(theta)]
%    V - Estimator variance (based on delta method approximation)
%    ESS - Effective Sample Size
%    eta1, eta2 - final values for optimal continuation probabilities
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
eta1 = 1;
eta2 = 1;
R_0 = 0;

dim_theta = length(p());
theta = zeros(dim_theta, M);

w = zeros(1,M);

I_k = false(1,M);
w_approx = zeros(1,M);
c_approx = zeros(1,M);
w_exact = zeros(1,M);
c_exact = zeros(1,M);

U = unifrnd(0,1,1,N);
update_factor = ceil(N/1000);

k = 0;
burning_in = true;
for i = 1:N
    % generate trial from the prior
    theta_trial = p();
    % simulate approximate data using these parameters
    start_t = toc;
    [D_s_approx, couple_arg_1, couple_arg_2, couple_arg_3] = s_approx(theta_trial);
    c_approx(i) = toc - start_t;
    theta(:,i) = theta_trial;
    % compute early accept/reject weight
    w_approx(i) = (rho_approx(D_s_approx) <= epsilon_approx);
    w(i) = w_approx(i);
    % compute continuation probility
    eta =  eta1*w(i) + eta2*(1- w_approx(i));
    % continue with probability eta
    if U(i) < eta
        % simulate exact model
        start_t = toc;
        D_s = s(theta_trial, couple_arg_1, couple_arg_2, couple_arg_3);
        c_exact(i) = toc - start_t;
        % update weights
        w_exact(i) = (rho(D_s) <= epsilon);
        w(i) = w(i) + (w_exact(i) - w(i))/eta;
        I_k(i) = true;
        if burning_in
            M = M-1;
            if M == 0
                burning_in = false;
            end
        end
    end
    
    if ~burning_in
        % recalculate values in 6.4--6.6
        rho_n = mean(w_approx(1:i));
        rho_k = mean(w_approx(I_k));
        % 6.4
        Ec = mean(c_approx(1:i));
        cp = (rho_n/rho_k)*mean(c_exact(I_k).*w_approx(I_k));
        cn = ((1-rho_n)/(1-rho_k))*mean(c_exact(I_k).*(1-w_approx(I_k)));
        % 6.5 
        ptp = (rho_n/rho_k)*mean(w_exact(I_k).*w_approx(I_k));
        pfp = (rho_n/rho_k)*mean((1-w_exact(I_k)).*w_approx(I_k));
        pfn = ((1-rho_n)/(1-rho_k))*mean(w_exact(I_k).*(1-w_approx(I_k)));
        % 6.6
        mu = sum(w(1:i).*f(theta(:,1:i)))/sum(w(1:i)) ;
        ptpf = (rho_n/rho_k)*mean(w_exact(I_k).*w_approx(I_k).*(f(theta(:,I_k))-mu).^2);
        pfpf = (rho_n/rho_k)*mean((1-w_exact(I_k)).*w_approx(I_k).*(f(theta(:,I_k))-mu).^2);
        pfnf = ((1-rho_n)/(1-rho_k))*mean(w_exact(I_k).*(1-w_approx(I_k)).*(f(theta(:,I_k))-mu).^2);
        % estimate optimal eta1 and eta2 (using f-depenent efficiency measure)
        R_0 = ptpf - pfpf;
        if R_0 > 0
            R_p = pfpf/(cp/Ec);
            R_n = pfnf/(cn/Ec);
            if max(R_p,R_n) <= R_0
                eta1 = sqrt(R_p/R_0);
                eta2 = sqrt(R_n/R_0);
            else
                eta1b = min(1,sqrt(R_p/R_0)/sqrt((1+pfnf/R_0)/(1+cn/Ec)));
                eta2b = min(1,sqrt(R_n/R_0)/sqrt((1+pfpf/R_0)/(1+cp/Ec)));
                % if phi(1,eta2) <= phi(eta1,1)
                if  (R_0 + pfpf +pfnf/eta2b)*(Ec + cp + eta2b*cn) <= (R_0 + pfpf/eta1b +pfnf)*(Ec + eta1b*cp + cn) 
                    eta1 = 1;
                    eta2 = eta2b;
                else
                    eta1 = eta1b;
                    eta2 = 1;
                end
            end
            % to avoid reducing to ABC rejection on the approximate model only.
            eta1 = max(eta1,1e-2);
            eta2 = max(eta2,1e-2);
        else
            eta1 = 1;
            eta2 = 1;
        end
    end
    
    if rem(i,update_factor)==0
        fprintf('%0.3f complete, eta1 = %0.3f ; eta2 = %0.3f ; R_0 = %0.3e \n', i/N, eta1, eta2, R_0);
    end
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

