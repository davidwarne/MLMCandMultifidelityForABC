%% Demonstration of Monte Carlo methods for
% approximate Bayesian computation 
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology


% initialise random number generator for reproducibility
rng(513,'twister');
% NOTE: parameters for the number of levels provided by external script
L = 6;
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
kmax = [4;30];
kmin = [1;10];

% functional of interest (i.e., computes mean of particular parameter)
f = @(x) x(1,:);


%% Set up ABC MLMC 
epsilon = zeros(L,1);
epsilon(1) = 1600;
epsilon(L) = 50;
% determine scale factor
m = exp(log(epsilon(1)/epsilon(L))/(L-1));
for l=(L-1):-1:2
    epsilon(l) = epsilon(l+1)*m;
end
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
tic;
% sequence of sample numbers
%N = [800;400;200;100;50]; % TODO: apply optimal choice
h = 1; % target RMSE 
N = ABCMLMCN(100,p,supp0,s,rho,epsilon,f,h)
C_mlmc_tune = toc;
% optimal N
%% Run and time ABC MLMC
fprintf('Running ABC MLMC...\n');
tic;
%[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,f,rho,epsilon)
[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,s,rho,epsilon,f);
C_mlmc = toc;
fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);

save(['MLMC_Rep_L',num2str(L),'_',num2str(m),'.mat']);
