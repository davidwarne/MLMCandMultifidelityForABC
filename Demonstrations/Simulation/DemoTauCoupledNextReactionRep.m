%% Demonstration of coupling Tau Leap with (Coupled) Next Reaction Method
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
for tau=0.01:0.01:0.1
rng(502,'twister');
h = figure;
% Build the repressilator  model
%[michment] = MichaelisMenten([0.001;0.005;0.01],100,100);
rep = Repressilator([1;1000;1;2;5;1],[0;0;0],[2;1;3]);
% simulate
[X_tau, t_tau, D, P] = TauLeapingMethod(rep,100,tau);
[X_ind, t_ind] = GillespieDirectMethod(rep,100);
[X_cpl, t_cpl] = CoupledNextReactionMethod(rep,100,D,P);

% Plot
Xn_tau = reshape([X_tau;X_tau],size(X_tau).*[1,2]); Xn_tau(:,end) = [];
tn_tau = reshape([t_tau;t_tau],[1,2*length(t_tau)]);tn_tau(1) = []; 

Xn_ind = reshape([X_ind;X_ind],size(X_ind).*[1,2]); Xn_ind(:,end) = [];
tn_ind = reshape([t_ind;t_ind],[1,2*length(t_ind)]);tn_ind(1) = []; 

Xn_cpl = reshape([X_cpl;X_cpl],size(X_cpl).*[1,2]); Xn_cpl(:,end) = [];
tn_cpl = reshape([t_cpl;t_cpl],[1,2*length(t_cpl)]);tn_cpl(1) = []; 

subplot(1,2,1);
plot(tn_cpl,Xn_cpl,tn_tau,Xn_tau,'k:','LineWidth',2);
xlim([0,100]); ylim([0,250]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
%legend({'S','E','C','P','tau-leap'});
title(sprintf('Coupled (tau %f)',tau));

subplot(1,2,2);
plot(tn_ind,Xn_ind,tn_tau,Xn_tau,'k:','LineWidth',2);
xlim([0,100]); ylim([0,250]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
legend({'M1','P1','M2','P2','M3','P3','tau-leap'});
title(sprintf('Uncoupled (tau %f)',tau))



end