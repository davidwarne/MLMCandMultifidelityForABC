%% Demonstration of Monte Carlo methods for approximate Bayesian computation 
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


%
% initialise random number generator for reproducibility
rng(513,'twister');

epsilon = 350;
tau = 0.04;
N = 10000;
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

%% Set up ABC Mulitfidelity
s_approx = @(k) GenerateApproxObservations(rep,[k_true(1:2);k;k_true(5:6)],X0,1,Obs_I,t,sig,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(rep,[k_true(1:2);k;k_true(5:6)],X0,1,Obs_I,t,sig,c1,c2,c3);
p = @() unifrnd(kmin,kmax);

%% Run and Time ABC Multifidelity
for i = 1:10
    fprintf('Running Adaptive ABC Multifidelity ...\n');
    tic;
    Ni = N*2^(i-1);
    Mi = Ni/10;
    [E_mf,V_mf,ESS_mf,Csim_mf,eta1,eta2,pairs] = ABCAdaptiveGradientMultifidelity(Ni,Mi,p,s_cpl,rho,epsilon,s_approx,rho,epsilon,f);
    C_mf = toc;
    fprintf('ABC Adaptive Multifidelity Completed in %f Sec\n',C_mf)
    save(['Bench_AdaptiveMF_Rep_epsilon',num2str(epsilon),'_tau',num2str(tau),'_',num2str(i),'.mat']);
end


