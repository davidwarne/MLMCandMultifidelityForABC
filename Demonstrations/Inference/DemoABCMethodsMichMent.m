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
f = @(x) x(1,:);
%% Set up ABC MLMC 
% discrepancy threshold sequence
epsilon = [40;20;10;5;2.5]; % TODO: vary as tuning param with L
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
% sequence of sample numbers
N = [800;400;200;100;50]; % TODO: apply optimal choice
%N = ABCMLMCN(100,p,supp0,f,rho,epsilon,0.5)

%% Run and time ABC MLMC
fprintf('Running ABC MLMC...\n');
tic;
%[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,f,rho,epsilon)
[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,s,rho,epsilon,f);
C_mlmc = toc;
fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);

%% Set up ABC Mulitfidelity
tau = 2;
s_approx = @(k) GenerateApproxObservations(michment,k,X0,1,[4],t,2,tau);
epsilon = 2.5;
p = @() unifrnd(kmin,kmax);
%N = 20000;
%tic;
%[eta1,eta2,p_fp,p_fn,p_tp,p_tn,ctilde,cp,cn,rho_dist,rho_dist_approx] = MultifidelityROC(N,p,f,rho,5,f_approx,rho,5);
%C_mfroc= toc;
%fprintf('ABC Multifidelity optimisation in %f Sec\n',C_mfroc)
% Prior sampler
N = 20000;
M = N/10;
%% Run and Time ABC Multifidelity
fprintf('Running Adaptive ABC Multifidelity ...\n');
tic;
[E_mf,V_mf,ESS_mf] = ABCAdaptiveMultifidelity(N,M,p,s,rho,epsilon,s_approx,rho,epsilon,f);
C_mf = toc;
fprintf('ABC Adaptive Multifidelity Completed in %f Sec\n',C_mf)

%% Set up ABC Rejection
% discrepancy threshold
epsilon = 2.5;
% Prior sampler
p = @() unifrnd(kmin,kmax);
% number of samples
N = 100;

%% Run and time ABC Rejection
fprintf('Running ABC Rejection...\n');
tic;
theta_rej = ABCRejectionSampler(N,p,s,rho,epsilon);
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
burnin = 100000;
thin = 1000;

%% Run and ABC MCMC
tic;
fprintf('Running ABC MCMC...\n');
theta0 = ABCRejectionSampler(1,p,f,rho,epsilon);
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
comp(:,4:6) = 1.95*sqrt(comp(:,4:6));

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
