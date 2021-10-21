%% Demonstration of Monte Carlo methods for
% approximate Bayesian computation 
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


%
% initialise random number generator for reproducibility
rng(513,'twister');
% NOTE: tau provided by external script
% tau = 0.01
% generate data from discrete sampling of a single realisation, 
% no observation error
k_true = [0.001;0.005;0.01]; 
X0 = [100;100;0;0];
t = [0;20;40;60;80];
[michment] = MichaelisMenten(k_true,X0(1),X0(2));
% assume no observation error for now.
Y_obs = GenerateObservations(michment,k_true,X0,1,[4],t,0);

% discrepancy function as a function of simulated data
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));

% Simulation as a function of k only
s = @(k) GenerateObservations(michment,k,X0,1,[4],t,2);
% prior support (uniform)
kmax = [0.003;0.015;0.05];
kmin = [eps;eps;eps];

% functional of interest (i.e., computes mean of particular parameter)
f = @(x) x(3,:);

%% Set up ABC Mulitfidelity
s_approx = @(k) GenerateApproxObservations(michment,k,X0,1,[4],t,2,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(michment,k,X0,1,[4],t,2,c1,c2,c3);
%epsilon = 200;
p = @() unifrnd(kmin,kmax);
%N = 20000;
%tic;
%[eta1,eta2,p_fp,p_fn,p_tp,p_tn,ctilde,cp,cn,rho_dist,rho_dist_approx] = MultifidelityROC(N,p,f,rho,5,f_approx,rho,5);
%C_mfroc= toc;
%fprintf('ABC Multifidelity optimisation in %f Sec\n',C_mfroc)
% Prior sampler
%N = 10000;
%M = N/10;
%% Run and Time ABC Multifidelity
fprintf('Running Adaptive ABC Multifidelity ...\n');
tic;
[E_mf,V_mf,ESS_mf,Csim_mf,eta1,eta2,pairs,~,~,ROC,Ec,cn,cp] = ABCAdaptiveGradientMultifidelity(N,M,p,s_cpl,rho,epsilon,s_approx,rho,epsilon,f);
C_mf = toc;
fprintf('ABC Adaptive Multifidelity Completed in %f Sec\n',C_mf)
save(['AdaptiveMF_MichMent_epsilon',num2str(epsilon),'_tau',num2str(tau),'N',num2str(N),'M',num2str(M),'smaller.mat']);

