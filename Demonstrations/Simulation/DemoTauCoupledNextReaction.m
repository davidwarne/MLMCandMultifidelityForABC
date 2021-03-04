%% Demonstration of coupling Tau Leap with (Coupled) Next Reaction Method
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% Build Michaelis-Menten model
[michment] = MichaelisMenten([0.001;0.005;0.01],100,100);
% simulate
[X_tau, t_tau, D, P] = TauLeapingMethod(michment,80,2);
[X_ind, t_ind] = GillespieDirectMethod(michment,80);
[X_cpl, t_cpl] = CoupledNextReactionMethod(michment,80,D,P);

% Plot
Xn_tau = reshape([X_tau;X_tau],size(X_tau).*[1,2]); Xn_tau(:,end) = [];
tn_tau = reshape([t_tau;t_tau],[1,2*length(t_tau)]);tn_tau(1) = []; 

Xn_ind = reshape([X_ind;X_ind],size(X_ind).*[1,2]); Xn_ind(:,end) = [];
tn_ind = reshape([t_ind;t_ind],[1,2*length(t_ind)]);tn_ind(1) = []; 

Xn_cpl = reshape([X_cpl;X_cpl],size(X_cpl).*[1,2]); Xn_cpl(:,end) = [];
tn_cpl = reshape([t_cpl;t_cpl],[1,2*length(t_cpl)]);tn_cpl(1) = []; 

subplot(1,2,1);
plot(tn_cpl,Xn_cpl,tn_tau,Xn_tau,'k:','LineWidth',2);
xlim([0,80]); ylim([0,100]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
%legend({'S','E','C','P','tau-leap'});
title('Coupled');

subplot(1,2,2);
plot(tn_ind,Xn_ind,tn_tau,Xn_tau,'k:','LineWidth',2);
xlim([0,80]); ylim([0,100]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
legend({'S','E','C','P','tau-leap'});
title('Uncoupled')