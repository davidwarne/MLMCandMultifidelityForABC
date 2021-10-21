function [bcrn] = TwoStepMAPKCascade(k,E0,X0,Y0,P10,P20)
%% Two step series MAPK enzymatic cascade
% Enzyme E phosphorates the substrate X to activate the protein and
% phosphotase P1 desphosphorylates the active protein. The activate X
% protein acts like an enzyme for the phosphoralation of substrate Y with
% phosphotase P2.
% 
% These reactions are:
%             k1
% R1:   X + E -> XE
%             k2
% R2:      XE -> X + E
%             k3
% R3:      XE -> X* + E
%             k4
% R4: X* + P1 -> X*P1
%             k5
% R5:    X*P1 -> X* + P1
%             k6
% R6:    X*P1 -> X + P1
%             k7
% R7: X* + Y -> X*Y
%             k8
% R8:    X*Y -> X* + Y
%             k9
% R9:    X*Y -> X* + Y*
%             k7
% R10: Y* + P2 -> Y*P2
%             k8
% R11:    Y*P2 -> Y* + P2
%             k9
% R12:    Y*P2 -> Y + P2
%    
% Inputs:
%    k - vector of kinetic rate parameters 
%    E0,X0,P10,Y0,P20 - initial copy numbers of enzyme, substrates and
%                       phosphotases
%   
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
% rate parameters 
bcrn.k = k;       
% initial copy numbers
bcrn.X0 = [E0;X0;0;0;P10;0;Y0;0;0;P20;0];
% number of reactions
bcrn.M = 12;               
% number of chemical species
bcrn.N = length(bcrn.X0);                         

%      1 2  3 4   5  6   7  8  9  10  11 
% X = [E,X,XE,X*,P1,X*P1,Y,YX*,Y*,P2,Y*P2]

% reactant stoichiometries
bcrn.nu_minus = [1,1,0,0,0,0,0,0,0,0,0;
                 0,0,1,0,0,0,0,0,0,0,0;
                 0,0,1,0,0,0,0,0,0,0,0;
                 0,0,0,1,1,0,0,0,0,0,0;
                 0,0,0,0,0,1,0,0,0,0,0;
                 0,0,0,0,0,1,0,0,0,0,0;
                 0,0,0,1,0,0,1,0,0,0,0;
                 0,0,0,0,0,0,0,1,0,0,0;
                 0,0,0,0,0,0,0,1,0,0,0;
                 0,0,0,0,0,0,0,0,1,1,0;
                 0,0,0,0,0,0,0,0,0,0,1;
                 0,0,0,0,0,0,0,0,0,0,1];
    
%      1 2  3 4   5  6   7  8  9  10  11 
% X = [E,X,XE,X*,P1,X*P1,Y,YX*,Y*,P2,Y*P2]
% product stoichiometries
bcrn.nu_plus = [0,0,1,0,0,0,0,0,0,0,0;
                1,1,0,0,0,0,0,0,0,0,0;
                1,0,0,1,0,0,0,0,0,0,0;
                0,0,0,0,0,1,0,0,0,0,0;
                0,0,0,1,1,0,0,0,0,0,0;
                0,1,0,0,1,0,0,0,0,0,0;
                0,0,0,0,0,0,0,1,0,0,0;
                0,0,0,1,0,0,1,0,0,0,0;
                0,0,0,1,0,0,0,0,1,0,0;
                0,0,0,0,0,0,0,0,0,0,1;
                0,0,0,0,0,0,0,0,1,1,0;
                0,0,0,0,0,0,1,0,0,1,0]; 
% stoichiometric matrix
bcrn.nu = bcrn.nu_plus - bcrn.nu_minus;
% propensity function
bcrn.a = @(X,k) k.*[X(1)*X(2);X(3);X(3);X(4)*X(5);X(6);X(6);X(4)*X(7);X(8);X(8);X(9)*X(10);X(11);X(11)];
