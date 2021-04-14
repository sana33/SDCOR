# LOF & LoOP

## Two state-of-the-art density-based techniques for outlier detection

In density-based methods, the local density of each object is calculated in a specific way and then is utilized to define the outlier scores. Given an object, the lower its local density compared to that of its neighbors, the more likely it is that the object is an outlier. Density around the points could be calculated by using many techniques, which most of them are distance-based. For example, Breunig et al. [1] propose a Local Outlier Factor (LOF) that uses the distance values of each object to its nearest neighbors to compute local densities. However, LOF has a drawback which is that the scores obtained through this approach are not globally comparable between all objects in the same dataset or even in different datasets. Kriegel et al. [2] introduce the Local Outlier Probability (LoOP), which is an enhanced version of LOF. LoOP gives each object a score in the interval [0,1], which is the probability of the object being an outlier and is widely interpretable among various situations.

[1] Breunig, Markus M., et al. "LOF: identifying density-based local outliers." Proceedings of the 2000 ACM SIGMOD international conference on Management of data. 2000.

[2] Kriegel, Hans-Peter, et al. "LoOP: local outlier probabilities." Proceedings of the 18th ACM conference on Information and knowledge management. 2009.

## Implementation details

The code combines the implementations of the two density-based techniques, as they both require a materialization matrix containing the _k_-nearest-neighbors and the following distances. You can follow the subsequent script with the suggested parameters as a template to use the `LOF_LoOP.m` function and obtain the LOF and LoOP results out of an arbitrary dataset:

```matlab
% Initializing input parameters
minPtsLB = 10; % LOF parameter
minPtsUB = 50; % LOF parameter
kStep= 2; % LOF parameter
kVal = 30; % LoOP parameter
lambda = 3; % LoOP parameter
OlNbCnd = 0; % Overlapping Neighborhood (OlNb) Condition for the materialization matrix
NSMethod = 0; % Neighbor Search Method for the knnsearch function [0:kdtree, 1:exhaustive] while OlNb is equal to 0
blkSzLim = 5e3; % memory size limit while OlNb equals 1 (our very own implementation)

% Mammography dataset
X = load('Mammography_(11183by6_260o).mat');
y = X.y; X = X.X;

[lofVals_Mammography,lofKmat_Mammography,ROC_LOF_Mammography,PR_LOF_Mammography,tElapsed_LOF_Mammography,LoOPvals_Mammography,ROC_LoOP_Mammography,...
    PR_LoOP_Mammography,tElapsed_LoOP_Mammography] = LOF_LoOP('Mammography',X,y,minPtsLB,minPtsUB,kStep,kVal,lambda,OlNbCnd,NSMethod,blkSzLim);
```


