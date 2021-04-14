# ORCA

## ORCA (Optimal Reciprocal Collision Avoidance)

Bay and Schwabacher [1] propose an optimized nested-loop algorithm based on the _k_ Nearest Neighbors (_k_-NN) distances among objects, that has a near-linear time complexity and is shortly named ORCA (Optimal Reciprocal Collision Avoidance). ORCA shuffles the input dataset in random order using a disk-based algorithm and processes it in blocks, as there is no need to load the entire data into the memory. It keeps looking for a set of the user-defined number of data points as potential anomalies, as well as for their anomaly scores. The cut-off value is set as the minimum outlierness score of the set, and it will get updates if there is a data point having a higher score in other blocks of data. For data points that obtain a lower score than the cut-off, they should be pruned; this pruning scheme will only expedite the process of distance computation in the case of data being ordered in an uncorrelated manner. ORCA's worst-case time-complexity is <img src="https://latex.codecogs.com/svg.image?O\left&space;(&space;n^2&space;\right&space;)" title="O\left ( n^2 \right )" />, and the I/O cost for the data accesses is quadratic. For the anomaly definition, it can use either the distance to the _k_-th nearest neighbor or the average distance of _k_-NN.

[1] Bay, Stephen D., and Mark Schwabacher. "Mining distance-based outliers in near linear time with randomization and a simple pruning rule." Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining. 2003.

## Implementation details

You can follow the subsequent script with the suggested parameters as a template to use the `EnLOF.m` function and obtain the required results out of an arbitrary dataset:

```matlab
% Mammography dataset
% "scores" contains the outliers indices provided by ORCA along with the subsequent outlier scores
% "labels" contains the outlier labels for all data elements; 0 for inliers, and 1 for outliers

scorTmp = zeros(numel(labels),1);
scorTmp(scores(:,1)) = scores(:,2);
scores = scorTmp;

[~,~,~,ROC_Mammography] = perfcurve(labels,scores,1);
[~,~,~,PR_Mammography] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('ORCA result for Mammography:\t\tROC = %0.3f\t\tPR = %0.3f\n\n',ROC_Mammography,PR_Mammography);
save('res_ORCA_Mammography.mat','ROC_Mammography','PR_Mammography');
```


