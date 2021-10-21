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

epsilon = 300;
tau = 0.5;
N = 10000;
% generate data from discrete sampling of a single realisation, 
_true = [0.001;0.001/120;0.18;0.001;0.001/22;0.3;0.0001;0.0001/110;0.2;0.001;0.001/22;0.3];
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
s_approx = @(k) GenerateApproxObservations(MAPK,[k_true(1);k(1);k(2);k_true(4);k(3);k(4);k_true(7);k(5);k(6);k_true(10);k(7);k(8)],X0,1,Obs_I,t,sig,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(MAPK,[k_true(1);k(1);k(2);k_true(4);k(3);k(4);k_true(7);k(5);k(6);k_true(10);k(7);k(8)],X0,1,Obs_I,t,sig,c1,c2,c3);

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
    save(['Bench_AdaptiveMF_MAPK_epsilon',num2str(epsilon),'_tau',num2str(tau),'_',num2str(i),'.mat']);
end


