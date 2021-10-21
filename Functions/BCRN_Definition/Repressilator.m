function [bcrn] = Repressilator(k,M0,P0)
%% Repressilator gene regulatory network
% Construction of a stochastic oscillator
% for m genes
%     alpha0 + alpha(K^n/[K^n+P_j^n])       
% G_i ------------------------------> M_i, 
%     beta                 beta       gamma
% M_i  ->  M_i + P_i,  P_i  -> 0, M_i  ->  0
%
% where j = (i+1 mod m) + 1
%
% Inputs:
%    k - vector of kinetic rate parameters and hill parameters
%    M0 - vector of initial populations of mRNA molecules  
%    P0 - vector of initial populations of proteins
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
% kinetic rate parameters and hill parameters k = [alpha0,alpha,K,n,beta,gamma]
bcrn.k = k;       
% number of reactions
bcrn.M = 4*length(M0); % for reactions per gene               
% number of chemical species
bcrn.N = 2*length(M0);                         
% reactant stoichiometries
base_nu_minus = [0,0;  % production of Mi
                 1,0;  % degradation of Mi
                 1,0;  % production of Pi
                 0,1]; % degradation of Pi 
bcrn.nu_minus = zeros(size(base_nu_minus)*length(M0));
for i=1:length(M0)
    bcrn.nu_minus((4*(i-1)+(1:4)),(2*(i-1)+(1:2))) = base_nu_minus;
end
    
% product stoichiometries
base_nu_plus = [1,0;
                0,0;
                1,1;
                0,0];
bcrn.nu_plus = zeros(size(base_nu_plus)*length(M0));
for i=1:length(M0)
    bcrn.nu_plus((4*(i-1)+(1:4)),(2*(i-1)+(1:2))) = base_nu_plus;
end
% stoichiometric matrix
bcrn.nu = bcrn.nu_plus - bcrn.nu_minus;
% initial copy numbers
bcrn.X0 = reshape([M0';P0'],[bcrn.N,1]);
% propensity function
bcrn.a = @(X,k) propensities(X,k);

function [a] = propensities(X,k)
%% caclulate propensities for arbitrary chain of genes
alpha0 = k(1);
alpha = k(2);
K = k(3);
n = k(4);
beta = k(5);
gamma = k(6);

N = length(X)/2;

a = zeros(N*4,1);
for i=1:N
    j = mod(i+N-2,N)+1;
    a(4*(i-1)+1) = alpha0 + alpha*(1.0/(1.0 + (X(2*(j-1)+2)/K)^n));
    a(4*(i-1)+2) = gamma*X(2*(i-1)+1);
    a(4*(i-1)+3) = beta*X(2*(i-1)+1);
    a(4*(i-1)+4) = beta*X(2*(i-1)+2);
end

