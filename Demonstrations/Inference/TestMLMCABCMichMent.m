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
L,epsilonL
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


%% Set up ABC MLMC 
epsilon = zeros(L,1);
epsilon(1) = 160;
epsilon(L) = epsilonL;
% determine scale factor
m = exp(log(epsilon(1)/epsilon(L))/(L-1));
for l=(L-1):-1:2
    epsilon(l) = epsilon(l+1)*m;
end
epsilon
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
tic;
% sequence of sample numbers
%N = [800;400;200;100;50]; % TODO: apply optimal choice
Ns = ABCMLMCN(100,p,supp0,s,rho,epsilon,f)
S = ceil(Ns(L))/Ns(L);
N = ceil(Ns*S)*16
C_mlmc_tune = toc;
% optimal N
%% Run and time ABC MLMC
fprintf('Running ABC MLMC...\n');
tic;
[E_mlmc,V_mlmc,F_mlmc,c_mlmc] = ABCMLMC(N,p,supp0,s,rho,epsilon,f);
C_mlmc = toc;
fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);
save(['MLMC_MichMent_epsilon',num2str(epsilonL),'_L',num2str(L),'smaller.mat']);
