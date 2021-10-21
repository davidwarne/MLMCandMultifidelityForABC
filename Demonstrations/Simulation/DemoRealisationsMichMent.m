%% Demonstration of multiple stochastic realisations
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
% initialise random number generator for reproducibility
rng(513,'twister');
h = figure;

% Build Michaelis-Menten model
k_true = [0.001;0.005;0.01]; 
X0 = [1000;1000;0;0];
t = [0;20;40;60;80];
[michment] = MichaelisMenten(k_true,X0(1),X0(2));
% assume no observation error for now.
Y_obs = GenerateObservations(michment,k_true,X0,1,[4],t,2);
rng(513,'twister');
% generate N realisations
N = 1000;
tic
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = TauLeapingMethod(michment,80,0.005);
end
toc/N
% plot samples with trasparant overlay (hint: for large N, make alpha smaller)
alpha = 0.05;
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    ha = plot(tn,Xn(1,:),'b','LineWidth',2); ha.Color(4) = alpha; 
    hb = plot(tn,Xn(2,:),'r','LineWidth',2); hb.Color(4) = alpha;
    hc = plot(tn,Xn(3,:),'Color',[237,177,32]/255,'LineWidth',2); hc.Color(4) = alpha; 
    hd = plot(tn,Xn(4,:),'Color',[127,47,142]/255,'LineWidth',2); hd.Color(4) = alpha;
end
%errorbar(t,Y_obs,2*ones(size(Y_obs)),'k.','LineWidth',2);
%errorbar(t,Y_obs-2,)
xlim([0,80]); ylim([0,1000]); legend({'$S_t$','$E_t$','$C_t$','$P_t$','$y_{obs}(t)$'});
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
box on