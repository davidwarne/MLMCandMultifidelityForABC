function [X,t] = GillespieDirectMethod(bcrn,T)
%% Gillespie Direct Method
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
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers

% block size (reallocate a new block every B updates)
B = 1000;
% initialise
X = zeros(length(bcrn.X0),B);
t = zeros(1,B);
X(:,1) = bcrn.X0;
t(1) = 0;

% count reactions
i = 1;
while true
    % compute propensities
    a = bcrn.a(X(:,i),bcrn.k);
    % sample exponential waiting time
    dt = exprnd(1/sum(a));
    % check if the simulation is finished
    if t(i) + dt <= T
        i = i + 1;
        % sample the next reaction event
        j = randsample(bcrn.M,1,true,a);
        
        % update copy numbers
        X(:,i) = X(:,i-1) + bcrn.nu(j,:)';
        t(i) = t(i-1) + dt;
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


