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
%epsilonL = 150;
%L = 3;
%tau = 0.5;
%N = 1000;
% generate data from discrete sampling of a single realisation, 
% no observation error 
k_true = [0.001;0.001/120;0.18;0.001;0.001/22;0.3;0.0001;0.0001/110;0.2;0.001;0.001/22;0.3];
X0  = [94;757; 0; 0;32;   0;567;  0; 0;32;   0];
t = linspace(0,200,100);
[MAPK] = TwoStepMAPKCascade(k_true,X0(1),X0(2),X0(7),X0(5),X0(10));

% assume only proteins are observable and additive error of sigma = 10
Obs_I = [4,9];
sig = 1;
Y_obs = GenerateObservations(MAPK,k_true,X0,1,Obs_I,t,sig);

% discrepancy function as a function of simulated data
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));

% Simulation as a function of k only

s = @(k) GenerateObservations(MAPK,[k_true(1);k(1);k(2);k_true(4);k(3);k(4);k_true(7);k(5);k(6);k_true(10);k(7);k(8)],X0,1,Obs_I,t,sig);
% prior support (uniform)
kmax = [k_true(1);1;k_true(4);1;k_true(7);1;k_true(10);1];
kmin = [0;0;0;0;0;0;0;0];

% functional of interest (i.e., computes mean of particular parameter)
f = @(x) x(7,:);

%% Set up ABC Mulitfidelity approximate model
% for the moment assume a fixed tau for all levels
sl_approx = @(k) GenerateApproxObservations(MAPK,[k_true(1);k(1);k(2);k_true(4);k(3);k(4);k_true(7);k(5);k(6);k_true(10);k(7);k(8)],X0,1,Obs_I,t,sig,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(MAPK,[k_true(1);k(1);k(2);k_true(4);k(3);k(4);k_true(7);k(5);k(6);k_true(10);k(7);k(8)],X0,1,Obs_I,t,sig,c1,c2,c3);
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
Ns = ABCAdaptiveGradientMultifidelityMLMCN(10000,1000,p,supp0,s_cpl,rho,epsilon,...
    s_approx,rho_approx,epsilon_approx,f);
%%
S = ceil(Ns(L))/Ns(L);
Ns = ceil(Ns*S)*1600
C_mfmlmc_tune = toc;

%% Run and time ABC MLMC
for i=1:10
    fprintf('Running ABC MF MLMC...\n');
    tic;
    Ni = Ns*2^(i-1);
    [E_mfmlmc,V_mfmlmc,F_mfmlmc] = ABCAdaptiveGradientMultifidelityMLMC(Ni,1000*ones(size(Ni)),p,supp0,s_cpl,rho,epsilon,...
    s_approx,rho_approx,epsilon_approx,f);
    C_mfmlmc = toc;
    fprintf('ABC MF MLMC Completed in %f Sec\n',C_mfmlmc);
    save(['Bench_MFMLMC_MAPK_epsilon',num2str(epsilonL),'_L',num2str(L),'_',num2str(i),'.mat']);
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

