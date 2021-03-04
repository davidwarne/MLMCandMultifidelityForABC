function [N]= ABCMLMCN(M,p,supp0,s,rho,epsilon,f,h)
%% Estimate optimal sample size sequence for MLMC based ABC
%
% Inputs:
%    M - the number of trial samples at each level determine optimal N
%    p - prior distribution sampler, 
%    supp0 - initial support region for sampling
%    s - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - a sequence of discrepancy acceptance thresholds
%    f - functional to estimate E[f(theta)] w.r.t. the ABC posterior measure
%    h - target root meansquare error
%
% Outputs:
%    E - Joint posterior mean estimate
%    V - estimator variances
%    F - Marginal CDF estimates (continuous functions)
%
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers


% initialise
L = length(epsilon);

theta = cell(L,1); 
supp = supp0;


C = 1e16;
k = length(supp.l);
v = zeros(L,1);
c = zeros(L,1);
N = zeros(L,1);
x = [];
for j=1:k
    x = [x;linspace(supp.l(j),supp.u(j),100)];
end
Fl = cell(k,1);
F = cell(k,1);
Finv = cell(k,1);

for l=1:L
    l
    %  ABC rejection step
    pl = @() p(supp0.l,supp0.u);
    supp_prev = supp;
    start_t = toc;
    theta{l} = ABCRejectionSampler(M,pl,s,rho,epsilon(l));
    % compute marginal eCDFs and support
    for j=1:k
        Fl{j} = @(t) ppval(interp1(x(j,:),ksdensity(theta{l}(j,:),x(j,:), ...
                   'Support','positive','BoundaryCorrection','reflection',...
                   'Function','cdf'),'pchip','pp'),t);
        supp.l(j) = min(theta{l}(j,:));
        supp.u(j) = max(theta{l}(j,:));
    end
    if l == 1
        % compute initial F and F^-1
        for j=1:k
           F{j} = @(t) ppval(interp1(x(j,:),Fl{j}(x(j,:)),'pchip','pp'),t); 
           [~,I] = unique(F{j}(x(j,:)));
           Finv{j} = @(u) ppval(interp1(F{j}(x(j,I)),x(j,I),'pchip','pp'),u); 
        end
        E = mean(f(theta{l}));
        V = (1/(M-1))*(mean(f(theta{l}).^2) - E.^2); 
        v(:,l) = V;
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
            Yl = @(t) Fl{j}(t) - Flm1(t);
            Fu = @(t) F{j}(t) + Yl(t);
            % monotonicity correction
            % F(s) = [u(s) +l(s)]/2
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
        % update estimators
        Pl = mean(f(theta{l}) - f(theta_lm1));
        E = E + Pl;
        v(l) = (1/(M-1))*(mean((f(theta{l}) - f(theta_lm1)).^2) - Pl.^2);
        V = V + v(l);
    end
    c(l) = toc - start_t;
end

for l=1:L
    N(l) = sqrt(v(l)/c(l))*sum((v.*c).^(1/2))/(h^2);
end
