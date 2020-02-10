This zip archive file contains two text file (in addition to this readme file).

1. Simulation Code.txt - This file contains R codes that implements the WL2Boost and PQLBoost algorithms. The main function that implements the WL2Boost 
algorithm is the boostingWL2 function. This function requires two other functions - fwd.lm1 and fwd.lme functions. The main function that implements the
PQLBoost algorithm is the boostingPQL function. This function also requires two other functions - fwd.lmer1 and estmixed functions. The annotation provided 
describes the role of this functions. The latter part of the file contain codes generating the dataset for the simulation.

2. Application Code.txt - This file contains codes implementing WL2Boost algorithms on real dataset. The functions implementing the WL2Boost algorithm are 
repeated as in Simulation Code.txt file. In addition, the file contains a 5-fold cross-validation function - boosting.CV.5(). The rest of the 
file contains data processing codes and calls to implement WL2Boost algorithms and necessary cross-validation.