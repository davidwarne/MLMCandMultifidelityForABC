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
rng(513,'twister');
h = figure;

% Build Two Step MAPK Cascade model
k_true = [0.001;0.001/120;0.18;0.001;0.001/22;0.3;0.0001;0.0001/110;0.2;0.001;0.001/22;0.3];
%      E   X   EX X* P1 X*P1 Y  YX*  Y* P2  Y*P2
X0  = [94;757; 0; 0;32;   0;567;  0; 0;32;   0];
t = linspace(0,200,50);
%                                   E     X    Y     P1    P2
[MAPK] = TwoStepMAPKCascade(k_true,X0(1),X0(2),X0(7),X0(5),X0(10));

% assume only proteins are observable and additive error of sigma = 10
Obs_I = [4,9];
sig = 10;
Y_obs = GenerateObservations(MAPK,k_true,X0,1,Obs_I,t,sig);
rng(513,'twister');
% generate N realisations
N = 1;
tic;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = GillespieDirectMethod(MAPK,200);
end
toc;
subplot(2,2,[1,2]);
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    plot(tn,Xn(1,:),'b','LineWidth',1); % E
    plot(tn,Xn(2,:),'r--','LineWidth',1);% X
    plot(tn,Xn(4,:),'r-','LineWidth',1);% X*
    plot(tn,Xn(5,:),'r:','LineWidth',1);% P_1
    plot(tn,Xn(7,:),'--','Color',[237,177,32]/255,'LineWidth',1); % Y
    plot(tn,Xn(9,:),'Color',[237,177,32]/255,'LineWidth',1); % Y*
    plot(tn,Xn(10,:),':','Color',[237,177,32]/255,'LineWidth',1); % P_2
    plot(tn,Xn(3,:),'--','Color',[127,47,142]/255,'LineWidth',1);% EX
    plot(tn,Xn(6,:),'r-.','LineWidth',1);% X*P1
    plot(tn,Xn(8,:),'--','Color',[253,141,60]/255,'LineWidth',1);% YX*
    plot(tn,Xn(11,:),'-.','Color',[237,177,32]/255,'LineWidth',1);% Y*P2
end
errorbar(t,Y_obs(1,:),sig*ones(size(Y_obs(1,:))),'k.','LineWidth',1);
errorbar(t,Y_obs(2,:),sig*ones(size(Y_obs(2,:))),'k.','LineWidth',1);
xlim([0,200]); ylim([0,800]); 
legend({'$E$','$X$','$X^*$','$P_1$','$Y$','$Y^*$','$P_2$','$XE$','$X^*P_1$','$YX^*$','$Y^*P_2$','$\mathbf{y}_{obs}(t)$'},'Location','northeastoutside');
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
box on

subplot(2,2,3);
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    plot(tn,Xn(1,:),'b','LineWidth',1); % E
    plot(tn,Xn(2,:),'r--','LineWidth',1);% X
    plot(tn,Xn(4,:),'r-','LineWidth',1);% X*
    plot(tn,Xn(5,:),'r:','LineWidth',1);% P_1
    plot(tn,Xn(7,:),'--','Color',[237,177,32]/255,'LineWidth',1); % Y
    plot(tn,Xn(9,:),'Color',[237,177,32]/255,'LineWidth',1); % Y*
    plot(tn,Xn(10,:),':','Color',[237,177,32]/255,'LineWidth',1); % P_2
    plot(tn,Xn(3,:),'--','Color',[127,47,142]/255,'LineWidth',1);% EX
    plot(tn,Xn(6,:),'r-.','LineWidth',1);% X*P1
    plot(tn,Xn(8,:),'--','Color',[253,141,60]/255,'LineWidth',1);% YX*
    plot(tn,Xn(11,:),'-.','Color',[237,177,32]/255,'LineWidth',1);% Y*P2
end
errorbar(t,Y_obs(1,:),sig*ones(size(Y_obs(1,:))),'k.','LineWidth',1);
errorbar(t,Y_obs(2,:),sig*ones(size(Y_obs(2,:))),'k.','LineWidth',1);
xlim([0,40]); ylim([0,200]); 
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
box on
subplot(2,2,4);
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    plot(tn,Xn(2,:),'r--','LineWidth',1);% X
    plot(tn,Xn(4,:),'r-','LineWidth',1);% X*
    plot(tn,Xn(5,:),'r:','LineWidth',1);% P_1
    plot(tn,Xn(7,:),'--','Color',[237,177,32]/255,'LineWidth',1); % Y
    plot(tn,Xn(9,:),'Color',[237,177,32]/255,'LineWidth',1); % Y*
    plot(tn,Xn(10,:),':','Color',[237,177,32]/255,'LineWidth',1); % P_2
    plot(tn,Xn(3,:),'--','Color',[127,47,142]/255,'LineWidth',1);% EX
    plot(tn,Xn(6,:),'r-.','LineWidth',1);% X*P1
    plot(tn,Xn(8,:),'--','Color',[253,141,60]/255,'LineWidth',1);% YX*
    plot(tn,Xn(11,:),'-.','Color',[237,177,32]/255,'LineWidth',1);% Y*P2
end
xlim([160,200]); ylim([0,200]); 
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
box on
print(h,'MAPKDat','-depsc','-painters')
