function [X,t] = CoupledNextReactionMethod(bcrn,T,D_obs,P_obs)
%% Coupled Next Reaction Method
% An exact stochastic simulation algorithm conditioned on observed
% coarse-grained Poisson processes
%
% Inputs:
%    bcrn  - a biochemical reaction network struct
%    T     - the end time of the simulation
%    D_obs - the coarse-grained (tau-leap) reaction-time interval lengths
%    P_obs - the Poisson number of reaction firings in each interval
%
%    Latter two inputs are outputs of:
%    [X_tau, t_tau, D_obs, P_obs] = TauLeapingMethod(bcrn, T, tau)
%
% Outputs:
%    X    -  time series of copy number vectors
%    t    -  vector of reaction times
%
% Authors:
%   Thomas P. Prescott[1] (prescott@maths.ox.ac.uk)
%   David J. Warne[2,3,4] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] Mathematical Institute, University of Oxford, UK
%   [2] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [3] Centre for Data Science, Queensland University of Technology, Autralia
%   [4] ARC Centre of Excellence for Mathematical and Statistical Frontiers


% block size (reallocate a new block every B updates)
B = 1000;
% initialise
X = zeros(length(bcrn.X0),B);
t = zeros(1,B);
X(:,1) = bcrn.X0;
t(1) = 0;
T_r = zeros(bcrn.M,1);

P = zeros(bcrn.M,1);
P_ctr = ones(bcrn.M,1);
P_lim = sum(P_obs,2)+1;
for m = 1:bcrn.M
    % couple
    conditional_poisson_process{m} = fine_grain_channel(D_obs(m,:), P_obs(m,:));
    % take the first of each unit-time exponential variate
    P(m) = conditional_poisson_process{m}(1);
end

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
        P_ctr(mu) = P_ctr(mu) + 1;
        if P_ctr(mu)<=P_lim(mu)
            P(mu) = conditional_poisson_process{mu}(P_ctr(mu)); % Take the next coupled exponential wait
        else
            P(mu) = P(mu) + exprnd(1); % Add an exponential variable if the coupled channel has ran out
        end
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
end


%% Coupling the coarse-grained Poisson process into an exact Poisson point process
function [firing_d] = fine_grain_interval(d_0, delta_d, num_observed_firings)
    firing_d = d_0 + delta_d*sort(rand(1, num_observed_firings));
end

function [T_r] = fine_grain_channel(D_r, P_r)
    
    D_r_cum = [0, cumsum(D_r)];
    P_r_cum = [0, cumsum(P_r)];

    % Initialise
    T_r = zeros(P_r_cum(end)+1, 1);    
    num_intervals = length(D_r);
   
    for i = 1:num_intervals
        % Start and length of reaction channel 'distance'
        d_0 = D_r_cum(i);
        delta_d = D_r(i);
        % Indices of the (unit rate) firings observed in this interval
        % (according to the coarse-grained Poission process)
        f_0 = P_r_cum(i);
        delta_f = P_r(i);
        
        % Place the firing points within the reaction channel interval
        T_r(f_0+1:f_0+delta_f) = fine_grain_interval(d_0, delta_d, delta_f);
    end
    
    % Ensure that the next firing is after the end of the coarse-grained
    % Poisson process (memoryless property)
    T_r(end) = D_r_cum(end) + exprnd(1);
end

