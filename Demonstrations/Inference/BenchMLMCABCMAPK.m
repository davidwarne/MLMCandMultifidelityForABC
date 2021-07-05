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
%L,epsilonL
%L = 6;
% generate data from discrete sampling of a single realisation, 
% no observation error (set up based on Prescott and Baker 2020)
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

%% Set up ABC MLMC 
epsilon = zeros(L,1);
epsilon(1) = 1600;
epsilon(L) = epsilonL;
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
Ns = ABCMLMCN(100,p,supp0,s,rho,epsilon,f)

S = ceil(Ns(L))/Ns(L);


Ns = ceil(Ns*S)*16;
C_mlmc_tune = toc;
% optimal N
%% Run and time ABC MLMC
for i=1:10
    fprintf('Running ABC MLMC...\n');
    tic;
    Ni = Ns*2^(i-1);
    [E_mlmc,V_mlmc,F_mlmc,c_mlmc] = ABCMLMC(Ni,p,supp0,s,rho,epsilon,f);
    C_mlmc = toc;
    fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);
    save(['Bench_MLMC_Rep_epsilon',num2str(epsilonL),'_L',num2str(L),'_',num2str(i),'.mat']);
end

