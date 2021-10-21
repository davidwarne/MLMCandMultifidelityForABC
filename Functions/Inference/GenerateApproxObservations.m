function [Y_obs, varargout] = GenerateApproxObservations(bcrn,k,X0,N,I,t,sigma,tau)
%% Generate Approximate Observations
% Generates N approximate Observations of state dimensions I of stochastic process X(t)
% with rate parameters k and initial condition X(0) = X0 = X0 using the tau-leap method. Observations are 
% taken at time points t(1),..,t(Nt) with observation error 
% Y_obs(t(i)) = X(t(i)) + zeta, where zeta ~ N(0,sigma^2)
%
% Inputs:
%    bcrn - a BCRN structure
%    k    - true kinetic rate parameters for the model
%    X0   - true initial condition of the BCRN
%    N    - number of sample paths in dataset
%    I    - vector of observable chemical species indices
%    t    - vector of observation times
%    sigma - standard deviation of observation error
%    tau   - discrete step in tau-leap method
% 
% Outputs:
%     a Table of observations Y_obs
%    optional outputs D and P representing the coupling from tau leap
%
%
% Authors:
%   Thomas P. Prescott[1] (tprescott@turing.ac.uk)
%   David J. Warne[2,3,4] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] The Alan Turing Institute, London, UK
%   [2] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [3] Centre for Data Science, Queensland University of Technology, Autralia
%   [4] ARC Centre of Excellence for Mathematical and Statistical Frontiers

Nt = length(t);
Y_obs = zeros(length(I),Nt,N);
bcrn.k = k;

save_coupling = false;
if nargout==4
    save_coupling = true;
    varargout{1} = [];
    varargout{2} = [];
elseif nargout==1
else
    warning('Request only one or four outputs from GenerateApproxObservations')
end

for i=1:N
    % stochastic simulation, saving coupling info if required
    if save_coupling
        [X_r,t_r,D_r,P_r] = TauLeapingMethod(bcrn,t(Nt),tau);
    else
        [X_r,t_r] = TauLeapingMethod(bcrn,t(Nt),tau);
    end
    
    % just store the state of selected species and times
    for j=1:Nt
        [J] = find(t_r <= t(j));
        Y_obs(:,j,i) = X_r(I,J(end));
    end
    % save coupling information if required
    if save_coupling
        varargout{1} = [varargout{1}, D_r];
        varargout{2} = [varargout{2}, P_r];
    end
end
% perturb with noise
zeta = normrnd(0,sigma,size(Y_obs));
Y_obs = Y_obs + zeta;

% Reshape coupled Poisson processes and save measurement noise for further
% coupling.
if save_coupling
    [m,d] = size(D_r);
    varargout{1} = reshape(varargout{1}, m, d, []);
    varargout{2} = reshape(varargout{2}, m, d, []);
    varargout{3} = zeta;
end
