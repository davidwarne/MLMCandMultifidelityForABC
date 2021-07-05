%% Demonstration of multiple stochastic realisations
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology


% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% get from PRL paper
theta = [2.37e-3;
    9.34e-2;
    7.57e-1;
    9.88e-1;
    3.40e-1;
    2.70e-1;
    1.11e0;
    1.30e-1;
    1.75e0;
    2.56e-1;
    4.51e0;
    8.21e-1];

EGF = 10;
NGF0 = 0;
fNGFR0 =01;
C3GI0 = 1;
MekI0 = 1;
ErkI0 = 100;

[EGFRsig] = EGFRSignalling(theta,EGF,NGF0,fNGFR0,C3GI0,ErkI0,MekI0);
%%
% generate N realisations
N = 1;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = ModifiedNextReactionMethod(EGFRsig,60);
end

% plot samples with trasparant overlay (hint: for large N, make alpha smaller)
alpha = 0.5;

%       1   2    3    4     5     6    7   8    9    10    11    12   13   14   15
% X = [EGF,NGF,fNGFR,bNGFR,bEGFR,RasA,P90,Raf1A,C3GI,C3GA,Rap1A,MekI,MekA,ErkI,ErkA]
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    ha = plot(tn,Xn([15],:)); %ha(1).Color(4) = alpha;ha(2).Color(4) = alpha;ha(3).Color(4) = alpha;ha(4).Color(4) = alpha;ha(5).Color(4) = alpha;
    %ha = plot(tn,Xn(15,:),'b','LineWidth',2); ha.Color(4) = alpha; 
   % hb = plot(tn,Xn(7,:),'r','LineWidth',2); hb.Color(4) = alpha;
    %hc = plot(tn,Xn(3,:),'Color',[237,177,32]/255,'LineWidth',2); hc.Color(4) = alpha; 
    %hd = plot(tn,Xn(4,:),'Color',[127,47,142]/255,'LineWidth',2); hd.Color(4) = alpha;
end
%xlim([0,80]); ylim([0,100]); 
legend({'Erk'});
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
