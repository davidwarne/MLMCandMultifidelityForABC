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

% create uniform joint prior
p = @() unifrnd(kmin,kmax);
tic;
% sequence of sample numbers
%% Run and time ABC Rejection
%for i=1:10
    fprintf('Running ABC Rejection...\n');
    tic;
    [theta_rej] = ABCRejectionSampler(Ni,p,s,rho,epsilon);
    E_rej = mean(f(theta_rej));
    V_rej = (1/(Ni-1))*(mean(f(theta_rej).^2) - E_rej.^2)
    C_rej = toc;
    fprintf('ABC Rejection Completed in %f Sec\n',C_rej);
    save(['Bench_Rej_MAPK_epsilon',num2str(epsilon),'_N',num2str(Ni),'.mat']);
%end

