%% Demonstration of multiple stochastic realisations
%
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% no observation error (set up based on Prescott and Baker 2020)
k_true = [1;1000;20;2;5;1]; % [alpha0,alpha,K,n,beta,gamma] 
X0 = [0;40;0;20;0;60];
t = [0;1;2;3;4;5;6;7;8;9;10];
[rep] = Repressilator(k_true,X0([1,3,5]),X0([2,4,6]));

% assume only proteins are observable and additive error of sigma = 10
Obs_I = [2,4,6];
sig = 10;
Y_obs = GenerateObservations(rep,k_true,X0,1,Obs_I,t,sig);
rng(502,'twister')
% generate N realisations
N = 1;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = GillespieDirectMethod(rep,10);
end

hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    plot(tn,Xn(1,:),'b--','LineWidth',2);  
    plot(tn,Xn(2,:),'b-','LineWidth',2); 
    plot(tn,Xn(3,:),'r--','LineWidth',2); 
    plot(tn,Xn(4,:),'r-','LineWidth',2);
    plot(tn,Xn(5,:),'--','Color',[237,177,32]/255,'LineWidth',2);
    plot(tn,Xn(6,:),'-','Color',[237,177,32]/255,'LineWidth',2);
end
errorbar(t,Y_obs(1,:),sig*ones(size(Y_obs(1,:))),'k.','LineWidth',2);
errorbar(t,Y_obs(2,:),sig*ones(size(Y_obs(2,:))),'k.','LineWidth',2);
errorbar(t,Y_obs(3,:),sig*ones(size(Y_obs(3,:))),'k.','LineWidth',2);
xlim([0,10]); ylim([0,400]); legend({'$M_{1,t}$','$P_{1,t}$','$M_{2,t}$','$P_{2,t}$','$M_{3,t}$','$P_{3,t}$','$\mathbf{y}_{obs}(t)$'},'Location','northeastoutside');
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
box on
