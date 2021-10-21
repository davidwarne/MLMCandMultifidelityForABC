function [Z,t,varargout] = TauLeapingMethod(bcrn,T,tau)
%% The Tau Leaping Method
% An approximate stochastic simulation algorithm
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    T    - the end time of the simulation
%    tau   - the timestep
% Outputs:
%    Z    -  time series of copy number vectors
%    t    -  vector of times    
%    vargout (optional) - save coarse-graining of reaction channel
%    unit-rate Poisson process
%  
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   Thomas P. Prescott[4] (tprescott@turing.ac.uk)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers
%   [4] The Alan Turing Institute, London, UK

% initialise
Nt = floor(T/tau);
Z = zeros(length(bcrn.X0),Nt+1);
t = zeros(1,Nt+1);
Z(:,1) = bcrn.X0;

save_coupling = false;
if nargout==4
    save_coupling = true;
elseif nargout==2
else
    warning('Request only two or four outputs from TauLeapingMethod')
end

if save_coupling
    varargout{1} = zeros(bcrn.M, Nt+1); % 'Distance' travelled for each reaction channel - speed a * timestep tau
    varargout{2} = zeros(bcrn.M, Nt+1); % Observations of (unit rate) firings over that distance
end

for i=1:Nt
    % compute propensities
    Z(Z(:,i) < 0,i) = 0;
    a = bcrn.a(Z(:,i),bcrn.k); a(a < 0 ) = 0;
    % generate poisson variates
    Y = poissrnd(a*tau);
    % update copy numbers
    Z(:,i+1) = Z(:,i) + (bcrn.nu') * Y;
    t(i+1) = t(i) + tau;
    
    if save_coupling
        varargout{1}(:,i) = a*tau; 
        varargout{2}(:,i) = Y;
    end
end
