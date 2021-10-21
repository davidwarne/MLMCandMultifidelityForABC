function [bcrn] = MichaelisMenten(k,S0,E0)
%% MichaelisMenten
% Construction of a Michaelis-Menten enzyme kinetic model
%       k(1)    k(2)        k(3)
% S + E -> C, C -> S + E, C -> E + P
%
% Inputs:
%    k - vector of kinetic rate parameters
%    S0 - initial population of substrate
%    E0 - initial population of enzymes
%
% Outputs:
%    a BCRN struct
%
% Author:
%   David J. Warne[1,2,3] (david.warne@qut.edu.au)
%   
% Affiliations:
%   [1] School of Mathematical Sciences, Queensland University of Technology, Autralia
%   [2] Centre for Data Science, Queensland University of Technology, Autralia
%   [3] ARC Centre of Excellence for Mathematical and Statistical Frontiers

bcrn = struct();
% kinetic rate parameters
bcrn.k = k;       
% number of reactions
bcrn.M = length(k);                
% number of chemical species
bcrn.N = 4;                        
% reactant stoichiometries
bcrn.nu_minus = [1,1,0,0;
                 0,0,1,0;
                 0,0,1,0];
% product stoichiometries
bcrn.nu_plus = [0,0,1,0;
                1,1,0,0;
                0,1,0,1];
% stoichiometric matrix
bcrn.nu = bcrn.nu_plus - bcrn.nu_minus;
% initial copy numbers
bcrn.X0 = [S0;E0;0;0];
% propensity function
bcrn.a = @(X,k) k.*[X(1)*X(2);X(3);X(3)];
