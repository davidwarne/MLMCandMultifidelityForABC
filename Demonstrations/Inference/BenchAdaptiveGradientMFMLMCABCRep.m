%% Demonstration of Monte Carlo methods for
% approximate Bayesian computation 
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


%
% initialise random number generator for reproducibility
rng(513,'twister');

% NOTE: epsilon, tau, N and M provided by external script
epsilonL = 350;
L = 3;
tau = 0.04;
%N = 1000;
% generate data from discrete sampling of a single realisation, 
% no observation error (set up based on Prescott and Baker 2020)
k_true = [1;1000;20;2;5;1]; % [alpha0,alpha,K,n,beta,gamma] 
X0 = [0;40;0;20;0;60];
t = [0;1;2;3;4;5;6;7;8;9;10];
[rep] = Repressilator(k_true,X0([1,3,5]),X0([2,4,6]));

% assume only proteins are observable and additive error of sigma = 10
Obs_I = [2,4,6];
sig = 1;
Y_obs = GenerateObservations(rep,k_true,X0,1,Obs_I,t,sig);

% discrepancy function as a function of simulated data
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));

% Simulation as a function of k only

s = @(k) GenerateObservations(rep,[k_true(1:2);k;k_true(5:6)],X0,1,Obs_I,t,sig);
% prior support (uniform)
kmax = [30;4];
kmin = [10;1];

% functional of interest (i.e., computes mean of particular parameter)
f = @(x) x(1,:);

%% Set up ABC Mulitfidelity approximate model
% for the moment assume a fixed tau for all levels
sl_approx = @(k) GenerateApproxObservations(rep,[k_true(1:2);k;k_true(5:6)],X0,1,Obs_I,t,sig,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(rep,[k_true(1:2);k;k_true(5:6)],X0,1,Obs_I,t,sig,c1,c2,c3);
%epsilon = 200;

%% Set up ABC MLMC epsilon sequence
epsilon = zeros(L,1);
epsilon(1) = 1600;
epsilon(L) = epsilonL;
% determine scale factor
m = exp(log(epsilon(1)/epsilon(L))/(L-1));
for l=(L-1):-1:2
    epsilon(l) = epsilon(l+1)*m;
end

% set up sequence of MF approximations
s_approx = cell(L,1);
rho_approx = cell(L,1);
epsilon_approx = epsilon;
for l=1:L
    s_approx{l} = sl_approx;
    rho_approx{l} = rho;
end
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
tic;
% estimate MLMC sequence
Ns = ABCAdaptiveGradientMultifidelityMLMCN(10000,100,p,supp0,s_cpl,rho,epsilon,...
    s_approx,rho_approx,epsilon_approx,f);
%%
S = ceil(Ns(L))/Ns(L);
Ns = ceil(Ns*S)*1600
C_mfmlmc_tune = toc;

%% Run and time ABC MLMC
for i=1:10
    fprintf('Running ABC MLMC...\n');
    tic;
    Ni = Ns*2^(i-1);
    [E_mfmlmc,V_mfmlmc,F_mfmlmc] = ABCAdaptiveGradientMultifidelityMLMC(Ni,100*ones(size(Ni)),p,supp0,s_cpl,rho,epsilon,...
    s_approx,rho_approx,epsilon_approx,f);
    C_mfmlmc = toc;
    fprintf('ABC MF MLMC Completed in %f Sec\n',C_mlmc);
    save(['Bench_MFMLMC_Rep_epsilon',num2str(epsilonL),'_L',num2str(L),'_',num2str(i),'.mat']);
end
%%% Run and Time ABC Multifidelity
%for i = 1:10
%    fprintf('Running Adaptive ABC Multifidelity ...\n');
%    tic;
%    Ni = N*2^(i-1);
%    Mi = Ni/10;
%    [E_mf,V_mf,ESS_mf,Csim_mf,eta1,eta2,pairs] = ABCAdaptiveGradientMultifidelity(Ni,Mi,p,s_cpl,rho,epsilon,s_approx,rho,epsilon,f);
%    C_mf = toc;
%    fprintf('ABC Adaptive Multifidelity Completed in %f Sec\n',C_mf)
%    save(['Bench_AdaptiveMF_Rep_epsilon',num2str(epsilon),'_tau',num2str(tau),'_',num2str(i),'.mat']);
%end
%

