# MATLAB Code implementations of multifidelity multilevel Monte Carlo for accelerated 
# approximate Bayesian inference for stochastic biochemical reaction networks

This repository contains reference MATLAB functions and scripts to demonstrate the multifiledity multilevel Monte Carlo approach for approximate Bayesian computation (MF-MLMC-ABC).

## Developers

David J. Warne$^{1,2,3}$ (david.warne@qut.edu.au), https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en

Thomas P. Prescott$^4$, https://scholar.google.com/citations?user=6fh76OQAAAAJ&hl=en

1. School of Mathematical Sciences, Faculty of Science, Queensland Univeristy of Technology, Australia
2. Centre for Data Science, Queensland University of Technology, Australia
3. ARC Centre of Excellence for Mathematical and Statistical Frontiers, Australia
4. The Alan Turing Institute, London, United Kingdom

## Citation Information

This code is provided as supplementary information to the paper,

David J Warne, Thomas P Prescott, Ruth E Baker, and Matthew J Simpson. Multifidelity multilevel Monte Carlo to accelerate approximate Bayesian parameter inference for partially observed stochastic processes. ArXiv preprint (TBA) 

## Licensing
This source code is licensed under the GNU General Public License Version 3.
Copyright (C) 2021 David J. Warne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Contents

```bash
The directory structure is as follows
|-- init.m                                  Adds all functions to the MATLAB Path
|-- Functions
    |-- BCRN_Definition                     Creation of BCRN data structures
        |-- MichaelisMenten.m
        |-- Repressilator.m
        |-- TwoStepMAPKCascade.m
    |-- Simulation                          Exact and approximate simulation schemes and coupling methods
        |-- GillespieDirectMethod.m
        |-- ModifiedNextReactionMethod.m
        |-- TauLeapingMethod.m
        |-- CorTauLeapingMethod.m
        |-- CoupledNextReactionMethod.m
    |-- Inference                           ABC methods based on Multifidelity and MLMC methods
        |-- GenerateObservations.m
        |-- GenerateApproximateObservations.m
        |-- GenerateCoupledObservations.m
        |-- ABCRejectionSampler.m
        |-- MultifidelityROC.m
        |-- ABCMultifidelity.m
        |-- ABCAdaptiveGradientMultifidelity.m
        |-- ABCMLMC.m
        |-- ABCMLMCN.m
        |-- ABCAdaptiveGradientMultifidelityMLMC.m
        |-- ABCAdaptiveGradientMultifidelityMLMCN.m
|-- Demostrations
    |-- Simulation                          Example realisation plots of models
        |-- DemoRealisationsMichMent.m
        |-- DemoRealisationsRep.m
        |-- DemoRealisationsMAPK.m
    |-- Inference                           Benchmark runs
        |-- BenchABCRejRep.m
        |-- BenchAdaptiveGradientMFABCRep.m
        |-- BenchMLMCABCRep.m
        |-- BenchAdaptiveGradientMFMLMCABCRep.m
        |-- BenchABCRejMAPK.m
        |-- BenchAdaptiveGradientMFABCMAPK.m
        |-- BenchMLMCABCMAPK.m
        |-- BenchAdaptiveGradientMFMLMCABCMAPK.m
```
## Usage

1. Start MATLAB
2. In MATLAB browse to the repository folder MLMCandMultifidelityForABC/
3. In the MATLAB command prompt, run 
   `>> init` 
   to set up the paths of all the example implementations.
4. Type in the name of any demo script in Demonstrations/Simulation/ 
   `>> DemoRealisationsMAPK` 
   generates Figure 7 in the paper.

5. Convergence data can be generated using scipts in Demonstrations/Inference/ 
   `>> BenchAdaptiveGradientMFMLMCABCRep.m`
   generates data for the MF-MLMC-ABC method on the repressilator model with $L= 5$, $\tau = 0.04$ and $\epsilon = 350$ as given in Figure 6B in the paper.

