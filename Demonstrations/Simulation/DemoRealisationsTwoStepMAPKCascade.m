%% Demonstration of multiple stochastic realisations
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology


% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% Build Two Step MAPK Cascade model
[MAPK] = TwoStepMAPKCascade([0.001;0.001/120;0.18;0.001;0.001/22;0.3;0.0001;0.0001/110;0.2;0.001;0.001/22;0.3],94,757,567,32,32);

% generate N realisations
N = 100;
tic;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    %[X_r{i},t_r{i}] = GillespieDirectMethod(MAPK,200);
    %[X_r{i},t_r{i}] = ModifiedNextReactionMethod(MAPK,200);
    [X_r{i},t_r{i}] = TauLeapingMethod(MAPK,200,1);
end
toc;
% plot samples with trasparant overlay (hint: for large N, make alpha smaller)
alpha = 0.05;
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    ha = plot(tn,Xn(1,:),'b','LineWidth',2); ha.Color(4) = alpha; 
    hb = plot(tn,Xn(4,:),'r','LineWidth',2); hb.Color(4) = alpha;
    hc = plot(tn,Xn(9,:),'Color',[237,177,32]/255,'LineWidth',2); hc.Color(4) = alpha; 
    %hd = plot(tn,Xn(4,:),'Color',[127,47,142]/255,'LineWidth',2); hd.Color(4) = alpha;
end
%xlim([0,80]); ylim([0,100]); 
legend({'E','X*','Y*'});
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
