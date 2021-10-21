%% Demonstration of Monte Carlo methods for approximate Bayesian computation 
%
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers


% initialise random number generator for reproducibility
rng(513,'twister');

Ni = 2048
epsilon = 350

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

% create uniform joint prior
p = @() unifrnd(kmin,kmax);
%% Run and time ABC Rejection
fprintf('Running ABC Rejection...\n');
tic;
[theta_rej] = ABCRejectionSampler(Ni,p,s,rho,epsilon);
E_rej = mean(f(theta_rej));
V_rej = (1/(Ni-1))*(mean(f(theta_rej).^2) - E_rej.^2)
C_rej = toc;
fprintf('ABC Rejection Completed in %f Sec\n',C_rej);
save(['Bench_Rej_Rep_epsilon',num2str(epsilon),'_N',num2str(Ni),'.mat']);

