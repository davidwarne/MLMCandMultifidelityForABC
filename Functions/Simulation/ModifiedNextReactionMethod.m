function [X,t] = ModifiedNextReactionMethod(bcrn,T)
%% Modified Next Reaction Method
% An exact stochastic simulation algorithm
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    T    - the end time of the simulation
% Outputs:
%    X    -  time series of copy number vectors
%    t    -  vector of reaction times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% block size (reallocate a new block every B updates)
B = 1000;
% initialise
X = zeros(length(bcrn.X0),B);
t = zeros(1,B);
X(:,1) = bcrn.X0;
t(1) = 0;
T_r = zeros(bcrn.M,1);
% generate M unit-time exponential variates
P = exprnd(1,[bcrn.M,1]);

% count reactions
i = 1;
while true
    % compute propensities
    a = bcrn.a(X(:,i),bcrn.k);
    % determine which reaction channel fires next
    dt = (P - T_r) ./ a;
    dt(a <= 0) = Inf;
    [delta,mu] = min(dt);
    if t(i) + delta <= T
        i = i + 1;
        %update copy numbers
        X(:,i) = X(:,i-1) + bcrn.nu(mu,:)';
        t(i) = t(i-1) + delta;
        T_r = T_r + a*delta;
        % update next reaction time for the firing channel
        P(mu) = P(mu) + exprnd(1);
        % allocate new block
        if mod(i,B) == 0
            X = [X,zeros(length(bcrn.X0),B)];
            t = [t,zeros(1,B)];
        end
    else
        % trim trailing zeros
        I = find(t == 0);
        t(I(2:end)) = [];
        X(:,I(2:end)) = [];
        return;
    end
end

