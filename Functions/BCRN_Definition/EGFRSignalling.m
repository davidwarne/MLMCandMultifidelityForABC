function [bcrn] = EGFRSignaling(k,EGF,NGF0,fNGFR0,C3GI0,ErkI0,MekI0)
%% Reduced model of Epidermal Growth Factor Receptor Signalling
%
% This is a stochastic version of the reduced model obtained by Transtrum and Qiu (2014)
% using Model Reduction by Manifold Boundaries. The full deterministic model includes
% 48 parameters and 15 ODEs, whereas the reduced model incldues 12 parameters
% and 8 ODEs. The stochatic varient that we utilise here includes 12 parameters,
% 15 chemical species and 14 reactions.

% These reactions are:
%                    k1
% R1:   NGF + bNGFR  -> fNGFR + bNGFR
%                    1
% R2:            EGF -> bEGFR
%                    k2
% R3:          bEGFR -> bEGFR + RasA
%                    k3
% R4:          bNGFR -> bNGFR + RasA
%                     1
% R5:     RasA + P90 -> P90
%                    k4
% R6:           RasA -> RasA + Raf1A
%                    k5/(Raf1A + k6)
% R7:          Raf1A -> 0
%                    k7
% R8:   bNGFR + C3GI -> bNGFR + C3GA
%                    k8
% R9:           C3GA -> C3GA + Rap1A
%                     1
% R10:  Raf1A + MekI -> Raf1A + MekA
%                    k9
% R11:         Rap1A -> Rap1A + MeKA
%                    k10
% R12:          ErkA -> 0
%                    k11
% R13:   MekA + ErkI -> MekA + ErkA
%                    k12
% R14:          ErkA -> ErKA + P90
%
% Inputs:
%    k - vector of effective parameters (as identified by the model reduction)
%    NGF0 - initial populations of Neuronal Growth Factor proteins
%    fNGFR0 - initial population of free Neuronal Growth Factor receptors
%    C3GI0, MeKI0, ErKI0 - initial populations of inactive, guanine nucleotide-releasing protein, 
%                          mitogen-activated protein kinase kinase, and extracellular signal-regulated kinase.
%
% Outputs:
%    a BCRN struct
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

bcrn = struct();
% rate parameters 
bcrn.k = k;       
% initial copy numbers
bcrn.X0 = [EGF;NGF0;fNGFR0;0;0;0;0;0;C3GI0;0;0;MekI0;0;ErkI0;0];
% number of reactions
bcrn.M = 14;               
% number of chemical species
bcrn.N = length(bcrn.X0);                         

%       1   2    3    4     5     6    7   8    9    10    11    12   13   14   15
% X = [EGF,NGF,fNGFR,bNGFR,bEGFR,RasA,P90,Raf1A,C3GI,C3GA,Rap1A,MekI,MekA,ErkI,ErkA]

% reactant stoichiometries
bcrn.nu_minus = [0,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
                 0,0,0,1,0,0,0,0,1,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
                 0,0,0,0,0,0,0,1,0,0,0,1,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                 0,0,0,0,0,0,0,0,0,0,0,0,1,1,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    
%       1   2    3    4     5     6    7   8    9    10    11    12   13   14   15
% X = [EGF,NGF,fNGFR,bNGFR,bEGFR,RasA,P90,Raf1A,C3GI,C3GA,Rap1A,MekI,MekA,ErkI,ErkA]
% product stoichiometries
bcrn.nu_plus =  [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,1,1,0,0,0,0,0,0,0,0,0;
                 0,0,0,1,0,1,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,1,0,1,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,1,0,0,0,0,1,1,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,1,1,0,0,0,0;
                 0,0,0,0,0,0,0,1,0,0,0,1,1,0,0;
                 0,0,0,0,0,0,0,0,0,0,1,0,1,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                 0,0,0,0,0,0,0,0,0,0,0,0,1,0,1;
                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,1]

% stoichiometric matrix
bcrn.nu = bcrn.nu_plus - bcrn.nu_minus;
% propensity function
bcrn.a = @(X,k) propensities(X,k);

function [a] = propensities(X,k)
%% caclulate propensities for reactions
a = zeros(14,1);

a(1) = k(1)*X(2)*X(3); % R1:   NGF + bNGFR  -> fNGFR + bNGFR
a(2) = X(1);           % R2:            EGF -> EFG + bEGFR
a(3) = k(2)*X(5);               % R3:          bEGFR -> bEGFR + RasA
a(4) = k(3)*X(4);               % R4:          bNGFR -> bNGFR + RasA
a(5) = X(6)*X(7);               % R5:     RasA + P90 -> P90
a(6) = k(4)*X(6);               % R6:           RasA -> RasA + Raf1A
a(7) = k(5)*X(8)/(X(8) + k(6)); % R7:          Raf1A -> 0
a(8) = k(7)*X(4)*X(9);          % R8:   bNGFR + C3GI -> bNGFR + C3GI + C3GA
a(9) = k(8)*X(10);              % R9:           C3GA -> C3GA + Rap1A
a(10) = X(8)*X(12);             % R10:  Raf1A + MekI -> Raf1A + MekI + MekA
a(11) = k(9)*X(11);             % R11:         Rap1A -> MeKA
a(12) = k(10)*X(15);            % R12:          ErkA -> 0
a(13) = k(11)*X(13)*X(14);      % R13:   MekA + ErkI -> MekA + ErkI + ErkA
a(14) = k(12)*X(15);            % R14:          ErkA -> ErKA + P90
