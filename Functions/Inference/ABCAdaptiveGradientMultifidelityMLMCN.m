function [Ns,c,v]= ABCAdaptiveGradientMultifidelityMLMCN(N,M,p,supp0,s_cpl,rho,epsilon,...
                          s_approx,rho_approx,epsilon_approx,f,h)
%%  Estimate the optimal sequence of sample sizes for Adaptive Multifidelity 
%   Multilevel Monte Carlo Sampler for approximate Bayesian computaion to compute
%   the posterior mean
%
% Inputs:
%    N - The number of trial weigthed sampeles at each level to determine optimal N
%    M - A number of burnin samples for continuation probability estimation
%    p - prior distribution sampler, 
%    supp0 - initial support region for sampling
%    s_cpl - function that generates simulated data give a parameters set
%            and coupling information (tau-leap reaction time-intevals, 
%            poission pocess firing numbers and obesration noise)
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - a sequence of discrepancy acceptance thresholds
%    s_approx - a sequence of functions that generate simulated approximate data 
%               give a parameters set
%    rho_approx - sequence of approximate discrepancy metric, treated as a 
%                 function of simulated data only
%    epsilon_approx - a sequence of pproximate discrepancy tresholds
%    f - functional to estimate E[f(theta)] w.r.t. the ABC posterior measure
%    h - if length(N) == 1, then h is the target RMSE.
%
% TODO: think about outputting optimal etas
%
% Outputs:
%    E - Joint posterior mean estimate
%    V - estimator variances
%    F - Marginal CDF estimates (continuous functions)
%
% Authors:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   Thomas P. Prescott[4] (prescott@maths.ox.ac.uk)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers
%   [4] Mathematical Institute, University of Oxford, UK

% initialise
L = length(epsilon)

theta = cell(L,1) 
weights = cell(L,1) 
v = zeros(L,1);
vabc = zeros(L,1);
c = zeros(L,1);
Ns = zeros(L,1);
supp = supp0
C = 1e16
k = length(supp.l)
x = []
for j=1:k
    x = [x;linspace(supp.l(j),supp.u(j),100)];
end
Fl = cell(k,1)
F = cell(k,1)
Finv = cell(k,1)

for l=1:L
    %  Multifidelity ABC sampling step
    pl = @() p(supp0.l,supp0.u);
     start_t = toc;
    [El,Vl,~,~,~,~,~,theta{l},weights{l}] = ABCAdaptiveGradientMultifidelity(N,M,pl,s_cpl,rho,epsilon(l),...
        s_approx{l},rho_approx{l},epsilon_approx(l),f);
    % compute marginal eCDFs and support
    for j=1:k
        j
        % TODO: determine the MF version of this (really the MF version of ksdensity)
        Fl{j} = @(t) ppval(interp1(x(j,:),ksdensity(theta{l}(j,:),x(j,:),'weights',weights{l}, ...
                   'Support','positive','BoundaryCorrection','reflection',...
                   'Function','cdf'),'pchip','pp'),t);       
        %lb = @(t) max(Fl(repmat(x(j,:)',[1,length(t)])) .*...
        %              (repmat(x(j,:)',[1,length(t)]) <= repmat(t,[length(x(j,:)),1]))...
        %                 - C*(repmat(x(j,:)',[1,length(t)]) > repmat(t,[length(x(j,:)),1])));
        %ub = @(t) min(Fl(repmat(x(j,:)',[1,length(t)])) .*...
        %              (repmat(x(j,:)',[1,length(t)]) >= repmat(t,[length(x(j,:)),1]))...
        %             + C*(repmat(x(j,:)',[1,length(t)]) < repmat(t,[length(x(j,:)),1])));
        %eFl{j} = @(t) ppval(interp1(x(j,:),lb(x(j,:))/2 + ub(x(j,:))/2,'pchip','pp'),t);
    end
    if l == 1
        % compute initial F and F^-1
        for j=1:k
            j
           %F{j} = @(t) ppval(interp1(x(j,:),eFl{j}(x(j,:)),'pchip','pp'),t); 
           F{j} = @(t) ppval(interp1(x(j,:),Fl{j}(x(j,:)),'pchip','pp'),t); 
           [~,I] = unique(F{j}(x(j,:)));
           Finv{j} = @(u) ppval(interp1(F{j}(x(j,I)),x(j,I),'pchip','pp'),u); 
        end
        % obtained directly from the MF sample at level 1
        E = El;
        V = Vl;
        v(l) = V;
    else
        % generate approximate coupled l-1 samples
        theta_lm1 = zeros(size(theta{l}));
        for j=1:k
            theta_lm1(j,:) = Finv{j}(Fl{j}(theta{l}(j,:)));
        end
        for j=1:k
            % compute marginal bias correction
            Flm1 = @(t) ppval(interp1(x(j,:),ksdensity(theta_lm1(j,:),x(j,:), ...
                       'Support','positive','BoundaryCorrection','reflection',...
                       'Function','cdf'),'pchip','pp'),t);
            Yl = @(t) Fl{j}(t) - Flm1(t); % maybe = @(t) eFl{j}(t) - Flm1(t)
            Fu = @(t) F{j}(t) + Yl(t);
            % monotonicity correction
            % F(x) = [u(x) +l(x)]/2
            % lb(t) = sup Fu((-infty,t])m ub(t) = inf Fu([t,infty))
            lb = @(t) max(Fu(repmat(x(j,:)',[1,length(t)])) .*...
                         (repmat(x(j,:)',[1,length(t)]) <= repmat(t,[length(x(j,:)),1]))...
                         - C*(repmat(x(j,:)',[1,length(t)]) > repmat(t,[length(x(j,:)),1])));
            ub = @(t) min(Fu(repmat(x(j,:)',[1,length(t)])) .*...
                         (repmat(x(j,:)',[1,length(t)]) >= repmat(t,[length(x(j,:)),1]))...
                         + C*(repmat(x(j,:)',[1,length(t)]) < repmat(t,[length(x(j,:)),1])));
            F{j} = @(t) ppval(interp1(x(j,:),lb(x(j,:))/2 + ub(x(j,:))/2,'pchip','pp'),t);
           [~,I] = unique(F{j}(x(j,:)));
            Finv{j} = @(u) ppval(interp1(F{j}(x(j,I)),x(j,I),'pchip','pp'),u);
        end
        % update estimators (need to use the weights also)
        B = f(theta{l}) - f(theta_lm1);
        sigma_w = sum(weights{l});
        sigma_wb = sum(weights{l}.*B);
        Pl = sigma_wb/sigma_w;
        E = E + Pl;
        cent_w = weights{l} - (sigma_w/N);
        cent_wb = weights{l}.*B - (sigma_wb/N);
        v(l)= sum((Pl*cent_w - cent_wb).^2)/sigma_w^2
        V = V + v(l);
    end
    c(l) = toc - start_t;
end

for l=1:L
    Ns(l) = sqrt(v(l)/c(l))*sum((v.*c).^(1/2));
end
