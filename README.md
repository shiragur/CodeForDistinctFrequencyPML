# Code for the experiments presented in the paper "Instance Based Approximations to Profile Maximum Likelihood"



Description of the files:

qsdpmlobj: Implements the objective function corresponding to the convex relaxation that approximates the PML.

ApproximatePMLACSS21: Code for implementing the heuristic version of Algorithm 2. Uses qsdpmlobj while solving the convex program.

PseudoPMLEntropyEstimation: Implements PseudoPML appraoch for entropy estimation using the algorithm in "ApproximatePMLACSS21" as subroutine to compute the PML part of the estimate.

PseudoPMLTest: Generates samples for different distributions and invokes the algorithm in file PseudoPMLEntropyEstimation.

