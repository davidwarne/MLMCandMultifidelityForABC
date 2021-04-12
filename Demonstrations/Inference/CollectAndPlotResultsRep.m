%% Demonstration of Monte Carlo methods for
% approximate Bayesian computation 
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');



%% results folder location
dir = 'matfiles_exp1/';
mlmc_results.L = []
mlmc_results.Ef = [];
mlmc_results.Var = []; 
mlmc_results.epsilon = [];
mlmc_results.compTime = [];

amf_results.tau = []
amf_results.Ef = [];
amf_results.Var = []; 
amf_results.epsilon = [];
amf_results.compTime = [];

NLmin_var = 10;
N_var = 10000;
M_var = 1000;
i = 1;
j = 1;
% collect results from HPC experiments
for epsilon_var = [500,450,350,300,250,200]
    for tau_var = [0.64,0.32,0.16,0.08,0.04,0.02,0.01,0.005]
        try
        [dir,'AdaptiveMF_Rep_epsilon',num2str(epsilon_var),'_tau',num2str(tau_var),'N',num2str(N_var),'M',num2str(M_var),'.mat']
        load([dir,'AdaptiveMF_Rep_epsilon',num2str(epsilon_var),'_tau',num2str(tau_var),'N',num2str(N_var),'M',num2str(M_var),'.mat']);
        amf_results.tau(j) = tau;
        amf_results.Ef(j) = E_mf;
        amf_results.Var(j) = V_mf; 
        amf_results.epsilon(j) = epsilon;
        amf_results.compTime(j) = C_mf;
        j = j+1;
        end
    end
    
    for L_var = [2,3,4,5,6,7,8,9]
        try
            load([dir,'MLMC_Rep_epsilon',num2str(epsilon_var),'_L',num2str(L_var),'NLmin',num2str(NLmin_var),'.mat']);
            mlmc_results.L(i) = L;
            mlmc_results.Ef(i) = E_mlmc;
            mlmc_results.Var(i) = V_mlmc; 
            mlmc_results.epsilon(i) = epsilonL;
            mlmc_results.compTime(i) = C_mlmc;
            i = i + 1;
        end
    end
end
%%
figure;
subplot(2,1,1);
plot(mlmc_results.L(mlmc_results.epsilon== 500),mlmc_results.compTime(mlmc_results.epsilon== 500),':+')
hold on;
for epsilon_var = [450,350,300,250,200]
plot(mlmc_results.L(mlmc_results.epsilon== epsilon_var),mlmc_results.compTime(mlmc_results.epsilon== epsilon_var),':+')   
end
ylim([0,10000]);

ylabel('$C$')
xlabel('$L$')
grid on
%set(gca,'XScale','Log');
legend({'$\epsilon = 500$','$\epsilon = 450$','$\epsilon = 350$','$\epsilon = 300$','$\epsilon = 250$','$\epsilon = 200$'})
subplot(2,1,2);
errorbar(mlmc_results.L(mlmc_results.epsilon== 500),mlmc_results.Ef(mlmc_results.epsilon== 500),sqrt(mlmc_results.Var(mlmc_results.epsilon== 500)),':+')
hold on;
for epsilon_var = [450,350,300,250,200]
errorbar(mlmc_results.L(mlmc_results.epsilon== epsilon_var),mlmc_results.Ef(mlmc_results.epsilon== epsilon_var),sqrt(mlmc_results.Var(mlmc_results.epsilon== epsilon_var)),':+')   
end
ylim([15,25]);
plot([2,9],[20,20],'--k');
ylabel('$\hat{f}\pm \sqrt{\mathrm{V}[\hat{f}]}$')
xlabel('$L$')
grid on;
%set(gca,'XScale','Log');

%% figure;
subplot(2,1,1);
plot(amf_results.tau(amf_results.epsilon== 500),amf_results.compTime(amf_results.epsilon== 500),':+')
hold on;
for epsilon_var = [450,350,300,250,200]
plot(amf_results.tau(amf_results.epsilon== epsilon_var),amf_results.compTime(amf_results.epsilon== epsilon_var),':+')
end
ylim([3000,8000]);
legend({'$\epsilon = 500$','$\epsilon = 450$','$\epsilon = 350$','$\epsilon = 300$','$\epsilon = 250$','$\epsilon = 200$'})
ylabel('$C$')
xlabel('$\tau$')
set(gca,'XScale','Log');
grid on
subplot(2,1,2);
errorbar(amf_results.tau(amf_results.epsilon== 500),amf_results.Ef(amf_results.epsilon== 500),sqrt(amf_results.Var(amf_results.epsilon== 500)),':+')
hold on;
for epsilon_var = [450,350,300,250,200]
errorbar(amf_results.tau(amf_results.epsilon== epsilon_var),amf_results.Ef(amf_results.epsilon== epsilon_var),sqrt(amf_results.Var(amf_results.epsilon== epsilon_var)),':+')   
end
ylim([16,22]);
plot([0.005,0.64],[20,20],'--k');
ylabel('$\hat{f}\pm \sqrt{\mathrm{V}[\hat{f}]}$')
xlabel('$\tau$')
grid on;
set(gca,'XScale','Log');
