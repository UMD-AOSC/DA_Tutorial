# DA_Tutorial
This open source project originated as a 'hands-on' tutorial for the RIKEN International School on Data Assimilation (RISDA2018).

The tutorial focuses on the Lorenz-63 (Lorenz 1963) model. The purpose is to assimilate observational data drawn from a 'true' driver
system to synchronize a response system (the data assimilation system). Some of the most popular data assimilation methods are made
available for experimentation, including:
Optimal Interpolation
3D-Var
Ensemble Kalman Filter
4D-Var
Particle Filter
Hybrid Filter

This is based on code developed for:
Penny, S.G., 2017: Mathematical foundations of hybrid data assimilation from a synchronization perspective
Chaos 27, 126801 (2017); https://doi.org/10.1063/1.5001819. https://aip.scitation.org/doi/full/10.1063/1.5001819

The second application currently available is the Modular Arbitrary Order Ocean Atmospehre Model (MAOOAM;
https://github.com/Climdyn/MAOOAM). This is an example of a coupled data assimilation system (CDA; the model is a 2-layer atmosphere
and 1-layer ocean), in which multiple scales must be addressed simultaneously. The code has been developed for and by the class of
AOSC658E at the University of Maryland to study CDA.

~Professor Stephen G. Penny
Currently: 
Research Scientist
NOAA Physical Sciences Laboratory (PSL) at the NOAA Earth Systems Research Laboratories (ESRL) in Boulder, CO, USA
and
Cooperative Institute for Research in Environmental Sciences at the University of Colorado Boulder

(Formerly: Professor, University of Maryland College Park)
