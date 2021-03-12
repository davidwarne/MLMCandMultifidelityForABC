%% Demonstration of Monte Carlo methods for
% approximate Bayesian computation 
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology


% initialise random number generator for reproducibility
rng(513,'twister');

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

s = @(k) GenerateObservations(rep,[k_true(1);k_true(2);k(1);k(2);k_true(5);k_true(6)],X0,1,Obs_I,t,sig);
% prior support (uniform)
kmax = [4;30];
kmin = [1;10];

% functional of interest (i.e., computes mean of particular parameter)
f = @(x) x(1,:);

%% Set up ABC Mulitfidelity
tau = 0.1;
s_approx = @(k) GenerateApproxObservations(rep,[k_true(1);k_true(2);k(1);k(2);k_true(5);k_true(6)],X0,1,Obs_I,t,sig,tau);
s_cpl = @(k,c1,c2,c3) GenerateCoupledObservations(rep,[k_true(1);k_true(2);k(1);k(2);k_true(5);k_true(6)],X0,1,Obs_I,t,sig,c1,c2,c3);
epsilon = 50;
p = @() unifrnd(kmin,kmax);
%N = 20000;
%tic;
%[eta1,eta2,p_fp,p_fn,p_tp,p_tn,ctilde,cp,cn,rho_dist,rho_dist_approx] = MultifidelityROC(N,p,f,rho,5,f_approx,rho,5);
%C_mfroc= toc;
%fprintf('ABC Multifidelity optimisation in %f Sec\n',C_mfroc)
% Prior sampler
N = 10000;
M = N/10;
%% Run and Time ABC Multifidelity
fprintf('Running Adaptive ABC Multifidelity ...\n');
tic;
[E_mf,V_mf,ESS_mf] = ABCAdaptiveMultifidelity(N,M,p,s_cpl,rho,epsilon,s_approx,rho,epsilon,f);
C_mf = toc;
fprintf('ABC Adaptive Multifidelity Completed in %f Sec\n',C_mf)

%% Set up ABC MLMC 
% discrepancy threshold sequence
L = 5; m = 2;
epsilon = zeros(L,1);
epsilon(L) = 50
for l=(L-1):-1:1
    epsilon(l) = epsilon(l+1)*m;
end
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
% sequence of sample numbers
%N = [800;400;200;100;50]; % TODO: apply optimal choice
 h = sqrt(V_mf); % target RMSE 
N = ABCMLMCN(100,p,supp0,s,rho,epsilon,f,h)
% optimal N
%% Run and time ABC MLMC
fprintf('Running ABC MLMC...\n');
tic;
%[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,f,rho,epsilon)
[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,s,rho,epsilon,f);
C_mlmc = toc;
fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);


%% Set up ABC Rejection
% discrepancy threshold
epsilon = 100;
% Prior sampler
p = @() unifrnd(kmin,kmax);
% number of samples
N = 100;

%% Run and time ABC Rejection
fprintf('Running ABC Rejection...\n');
tic;
%theta_rej = ABCRejectionSampler(N,p,s,rho,epsilon);
[theta_rej,r] = ABCRejectionSamplerQuant(N,p,s,rho,0.05);

E_rej = mean(f(theta_rej));
V_rej = (1/(N-1))*(mean(f(theta_rej).^2) - E_rej.^2);
C_rej = toc;
fprintf('ABC Rejection Completed in %f Sec\n',C_rej);


%% Set up ABC MCMC 
% discrepancy threshold
epsilon = 2.5;
% Prior sampler and PDF
p = @() unifrnd(kmin,kmax);
p_pdf = @(k) prod(unifpdf(k,kmin,kmax));
% create proposal Kernel sampler and PDF
Sigma = diag((0.05*(kmax-kmin)).^2);
K = @(k) mvnrnd(k,Sigma)';
K_pdf = @(k_n,k_p) mvnpdf(k_n,k_p,Sigma);
% number of steps in the Markov Chain
T = 500000;
burnin = 1000;
thin = 1;%000;

%% Run and ABC MCMC
tic;
fprintf('Running ABC MCMC...\n');
theta0 = ABCRejectionSampler(1,p,s,rho,epsilon);
theta_mcmc = ABCMCMCSampler(T,p_pdf,K,K_pdf,s,rho,epsilon,theta0);
E_mcmc = mean(f(theta_mcmc(:,burnin:thin:T)));
V_mcmc = (1/(N-1))*(mean(f(theta_mcmc(:,burnin:thin:T)).^2) - E_mcmc.^2);
C_mcmc = toc;
fprintf('ABC MCMC Completed in %f Sec\n',C_mcmc);

%% Set up ABC SMC
% discrepancy threshold sequence
epsilon = [40;20;10;5;2.5];
% Prior sampler and PDF
p = @() unifrnd(kmin,kmax);
p_pdf = @(k) prod(unifpdf(k,kmin,kmax));
% create proposal Kernel sampler and PDF
Sigma = diag((0.05*(kmax-kmin)).^2);
K = @(k) mvnrnd(k,Sigma)';
K_pdf = @(k_n,k_p) mvnpdf(k_n,k_p,Sigma);
% number of particles
N = 100;

%% Run and time ABC SMC
fprintf('Running ABC SMC...\n');
tic;
[theta_smc,W] = ABCSMCSampler(N,p,p_pdf,K,K_pdf,s,rho,epsilon);
E_smc = mean(f(theta_smc(:,:,end)),2);
V_smc = (1/(N-1))*(mean(f(theta_smc(:,:,end)).^2,2) - E_smc.^2);
C_smc = toc;
fprintf('ABC SMC Completed in %f Sec\n',C_smc);

%% collect results
comp = [E_mlmc,V_mlmc,C_mlmc;
        E_mf,V_mf,C_mf;
        E_rej,V_rej,C_rej;
        E_mcmc,V_mcmc,C_mcmc;
        E_smc,V_smc,C_smc]
comp(:,2) = 1.95*sqrt(comp(:,2));

%% plot Markov Chain transient behaviour
for i=1:3
    figure;
    plot(theta_mcmc(i,:),'b','LineWidth',2);xlabel('t');ylabel(['k_',num2str(i)])
    hold on;
    plot([0,length(theta_mcmc(i,:))],[k_true(i),k_true(i)],'--k','LineWidth',2)
end

% plot marginal densities 
ks = [linspace(0,kmax(1),1000);linspace(0,kmax(2),1000);linspace(0,kmax(3),1000)];
for i=1:3
     figure;
     hold on;        
     plot(ks(i,2:end),diff(F_mlmc{i}(ks(i,:)))./diff(ks(i,:)));
     ksdensity(theta_rej(i,:),'Support','positive','BoundaryCorrection','reflection'); 
     ksdensity(theta_mcmc(i,burnin:thin:T),'Support','positive','BoundaryCorrection','reflection');
     ksdensity(theta_smc(i,:,end),'Support','positive','BoundaryCorrection','reflection'); 
     ylabel(['p_\epsilon (k_',num2str(i),')']);xlabel(['k_',num2str(i)])
     xlim([ks(i,1),ks(i,end)]);
end