function [E,V,ESS,c_sim,eta1,eta2,pairs] = ABCAdaptiveGradientMultifidelity(N,M,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx,varargin)
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
%    pairs - number of pairs computed
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

%-----------------------------------%
% Start with burn-in
[theta_burnin, w_approx, w_exact, c_approx, c_exact] = generate_burnin(M,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx);

% Keep track of simulation burden
c_sim = sum(c_approx) + sum(c_exact);

% Get post burn-in stats for an initial guess of optimal eta1, eta2

% Expected approx sim cost
Ec = mean(c_approx);
% Probability that approx is positive
p_pos = mean(w_approx);
% Probability that approx is positive GIVEN it's been continued
p_pos_given_cont = p_pos;

% Mean exact sim costs, split across positives and negatives
cp = mean(c_exact .* w_approx);
cn = mean(c_exact .* (1-w_approx));

% Save acceptances (with which we can calculate mu)
% Use w_exact because eta=1 means we ignore w_approx for acceptances
theta = theta_burnin(:,w_exact~=0);
w = w_exact(w_exact~=0);
Z = sum(w);

% Estimate functional mean, if required
try
    f = varargin{1};
catch
    f = @(theta) ones(1, size(theta,2));
end
mu = sum(w.*f(theta))/Z;

% Receiver Operator Characteristics
ROC_MAT = get_burnin_ROC(theta_burnin, w_approx, w_exact, varargin{:});

% Start stepping optimal eta
delta_evolve = 10;
[eta1, eta2] = evolve_eta(1, 1, p_pos, p_pos_given_cont, mu, ROC_MAT, Ec, cp, cn, delta_evolve);

%-----------------------------------%
% Begin main part of algorithm

U = unifrnd(0,1,1,N-M);
update_factor = ceil(N/100);

k=0;
for i = 1:(N-M)
    
    % generate trial from the prior
    theta_trial = p();
    
    % simulate approximate data using these parameters
    start_t = toc;
    [D_s_approx, couple_arg_1, couple_arg_2, couple_arg_3] = s_approx(theta_trial);
    c_approx_trial = toc - start_t;
    c_sim = c_sim + c_approx_trial;
        
    % update Ec (average approximate simulation time)
    Ec = Ec + (c_approx_trial - Ec)/(M+i);
        
    % compute early accept/reject weight
    w_approx_trial = (rho_approx(D_s_approx) <= epsilon_approx);
    w_trial = w_approx_trial;
    
    % update rho
    p_pos = p_pos + (w_approx_trial - p_pos)/(M+i);

    % compute continuation probility
    eta =  eta1*w_approx_trial + eta2*(1- w_approx_trial);
    
    % continue with probability eta
    if U(i) < eta
        k = k+1;
        
        % update Probability that approx is positive GIVEN it's been continued
        p_pos_given_cont = p_pos_given_cont + (w_approx_trial - p_pos_given_cont)/(M+k);
        
        % simulate exact model
        start_t = toc;
        D_s = s(theta_trial, couple_arg_1, couple_arg_2, couple_arg_3);
        c_exact_trial = toc - start_t;
        c_sim = c_sim + c_exact_trial;
        
        % update cp and cn
        cp = cp + (c_exact_trial*w_approx_trial - cp)/(M+k);
        cn = cn + (c_exact_trial*(1-w_approx_trial) - cn)/(M+k);
        
        % update weights
        w_exact_trial = (rho(D_s) <= epsilon);
        w_trial = w_trial + (w_exact_trial - w_trial)/eta;
    end
    
    % Store theta, w if non-zero weight
    % Also update mean of functional if changed
    if w_trial ~= 0
        theta = [theta, theta_trial];
        w = [w, w_trial];
        Z = Z + w_trial;
        mu = mu + (f(theta_trial) - mu)*(w_trial/Z);
    end
    
    % Update ROC_MAT
    if U(i) < eta
        dROC = get_burnin_ROC(theta_trial, w_approx_trial, w_exact_trial, varargin{:});
        ROC_MAT = ROC_MAT + (dROC - ROC_MAT)/(M+k);
    end
    
    [eta1, eta2] = evolve_eta(eta1, eta2, p_pos, p_pos_given_cont, mu, ROC_MAT, Ec, cp, cn, delta_evolve);
    
    if rem(i,update_factor)==0
        fprintf('%0.3f complete, eta1 = %0.3f ; eta2 = %0.3f ; Z = %0.3f ; pairs = %d \n', (M+i)/N, eta1, eta2, Z/(M+i), M+k);
    end
end

pairs = M + k;

F = f(theta);

% compute Multifidelity estimator
Sigma_w = sum(w);
Sigma_wf = sum(w.*F);
E = Sigma_wf/Sigma_w;

% compute Variance approximation
centred_w = w - (Sigma_w/N);
centred_wf = w.*F - (Sigma_wf/N);

V = sum((E*centred_w - centred_wf).^2)/Sigma_w^2;

% effective sample size (proportional to ESS \propto 1/V)
ESS = (Sigma_w^2)/sum(w.*w);


end

function [theta, w_approx, w_exact, c_approx, c_exact] = generate_burnin(M,p,s,rho,epsilon,s_approx,rho_approx,epsilon_approx)
%% Burn-in 

    % initialise
    dim_theta = length(p());
    theta = zeros(dim_theta, M);
    w_approx = zeros(1,M);
    c_approx = zeros(1,M);
    w_exact = zeros(1,M);
    c_exact = zeros(1,M);

    for i = 1:M
        % generate trial from the prior
        theta_trial = p();
        % simulate approximate data using these parameters
        start_t = toc;
        [D_s_approx, couple_arg_1, couple_arg_2, couple_arg_3] = s_approx(theta_trial);
        c_approx(i) = toc - start_t;
        theta(:,i) = theta_trial;
        % compute early accept/reject weight
        w_approx(i) = (rho_approx(D_s_approx) <= epsilon_approx);
        % simulate exact model
        start_t = toc;
        D_s = s(theta_trial, couple_arg_1, couple_arg_2, couple_arg_3);
        c_exact(i) = toc - start_t;
        % update weights
        w_exact(i) = (rho(D_s) <= epsilon);    
    end
end

function [ROC_MAT] = get_burnin_ROC(theta, w_approx, w_exact, varargin)
    %% ROC MATRIX
    % Produce a matrix to estimate:
    % ROC_MAT =
    % [sum(1 & TP) sum(1 & FP) sum(1 & FN) sum(1 & TN);
    %  sum(F & TP) sum(F & FP) sum(F & FN) sum(F & TN);
    %  sum(F^2 & TP) sum(F^2 & FP) sum(F^2 & FN) sum(F^2 & TN)];
    %
    % Designed so that [mu^2 -2*mu 1] * ROC_MAT = 
    % [sum((F-mu)^2 & TP) sum((F-mu)^2 & FP) sum((F-mu)^2 & FN) sum((F-mu)^2 & TN)]
    %
    % If no function supplied, ROC_MAT fixed such that
    % [mu^2 -2*mu 1] * ROC_MAT = [sum(TP) sum(FP) sum(FN) sum(TN)]
    
    M = length(w_exact);
    mf_mat = [w_approx.*w_exact; w_approx.*(1-w_exact); (1-w_approx).*w_exact; (1-w_approx).*(1-w_exact)];
    
    if nargin==4
        f = varargin{1};
        moment_mat = [ones(1,M); f(theta); f(theta).^2];
    elseif nargin==3
        moment_mat = [zeros(2,M); ones(1,M)];
    else
        error("get_ROC doesn't make sense without 3 or 4 inputs")
    end
    ROC_MAT = (moment_mat * mf_mat')/M;
end

function [eta1, eta2] = evolve_eta(eta1, eta2, rho_n, rho_k, mu, ROC_MAT, Ec, cp, cn, delta)
    mu_v = [mu^2 -2*mu 1];
    ptp = (rho_n/rho_k) * (mu_v * ROC_MAT(:,1));
    pfp = (rho_n/rho_k) * (mu_v * ROC_MAT(:,2));
    pfn = ((1-rho_n)/(1-rho_k)) * (mu_v * ROC_MAT(:,3));
    
    cp = (rho_n/rho_k) * cp;
    cn = ((1-rho_n)/(1-rho_k)) * cn;
    
    % estimate optimal eta1 and eta2 (using f-dependent efficiency measure)
    R0 = ptp - pfp;
    grad1 = pfp*cn*eta2/eta1 - pfn*cp*eta1/eta2 - cp*R0*eta1 + pfp*Ec/eta1;
    grad2 = pfn*cp*eta1/eta2 - pfp*cn*eta2/eta1 - cn*R0*eta2 + pfn*Ec/eta2;
    
    grad_scale =  (Ec + cp + cn)*(mu^2);    
    eta1 = min(1, eta1 * exp(delta*grad1/grad_scale));
    eta2 = min(1, eta2 * exp(delta*grad2/grad_scale));
end
