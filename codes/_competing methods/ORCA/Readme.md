# ORCA

## ORCA (Optimal Reciprocal Collision Avoidance)

Bay and Schwabacher [1] propose an optimized nested-loop algorithm based on the _k_ Nearest Neighbors (_k_-NN) distances among objects, that has a near-linear time complexity and is shortly named ORCA (Optimal Reciprocal Collision Avoidance). ORCA shuffles the input dataset in random order using a disk-based algorithm and processes it in blocks, as there is no need to load the entire data into the memory. It keeps looking for a set of the user-defined number of data points as potential anomalies, as well as for their anomaly scores. The cut-off value is set as the minimum outlierness score of the set, and it will get updates if there is a data point having a higher score in other blocks of data. For data points that obtain a lower score than the cut-off, they should be pruned; this pruning scheme will only expedite the process of distance computation in the case of data being ordered in an uncorrelated manner. ORCA's worst-case time-complexity is <img src="https://latex.codecogs.com/svg.image?O\left&space;(&space;n^2&space;\right&space;)" title="O\left ( n^2 \right )" />, and the I/O cost for the data accesses is quadratic. For the anomaly definition, it can use either the distance to the _k_-th nearest neighbor or the average distance of _k_-NN.

In ORCA, the parameter _k_ denotes the number of nearest neighbors, which by using higher values for that, the execution time will also increase. Here, the suggested value for _k_, equal to 5, is utilized in all of our experiments. The parameter _N_ specifies the maximum number of anomalies to be reported. If _N_ is set to a small value, then ORCA increases the running cut-off quickly, and therefore, more searches will be pruned off, which will result in a much faster runtime. Hence, as the correct number of anomalies is not assumed to be foreknown in the algorithm, we set <img src="https://latex.codecogs.com/svg.image?N=\frac{n}{8}" title="N=\frac{n}{8}" />, as a reasonable value, where _n_ stands for the cardinality of the input data.

However, ORCA does not report an anomaly score for the rest of the data, and this is as long as for computing the AUC values, there is an indispensable need to report an anomaly score for every instance; thus, we set a score equal to zero for other non-anomaly reported objects.

[1] Bay, Stephen D., and Mark Schwabacher. "Mining distance-based outliers in near linear time with randomization and a simple pruning rule." Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining. 2003.

## Implementation details

**The implementation is in C++ and provided by the genuine authors through [this link](http://www.stephenbay.net/orca/).**

In the case of datasets in MAT format, you shall convert them to a format acceptable by ORCA. For this matter, first, you should save them through MATLAB in ASCII format; this can be done using the `dlmwrite()` or the `writematrix()` functions. Then you must employ `dprep.exe` given by the ORCA authors along with the mentioned necessities to change the ASCII data into a binary file required by ORCA.

Finally, for running the code, it will be only required to use `orca.exe` with the suggested parameters to obtain the anomaly scores. You can utilize the `ptime.exe` executable code to calculate the runtime for each command and save the command output in a specific file by adding `> orca_output.txt` to the end of the command line; outlier scores and the execution time can be elicited out of this file.

You can follow the subsequent script with the suggested parameters as a template to use the ORCA implementation code and obtain the required results out of an arbitrary dataset:

```matlab
%%% Mammography dataset
```

```matlab
%% converting the MAT dataset into a binary format acceptable by ORCA

1> in MATLAB:

load('Mammography_(11183by6_260o).mat');
dlmwrite('Mammography',X,'precision','%.15f'); dlmwrite('labels',y);

2> create the "Mammography.fields" file with the following content:

attrib01: continuous.
attrib02: continuous.
attrib03: continuous.
attrib04: continuous.
attrib05: continuous.
attrib06: continuous.

3> in command window:

dprep.exe Mammography Mammography.fields Mammography.bin -rand -snone -cleanf
```

```matlab
%% running the ORCA C++ executable code in command window

C:\ptime.exe orca.exe Mammography.bin Mammography.bin Mammography.weights -n 1397 > Mammography_ORCA.comOut
```

```matlab
%% gaining the AUC outcomes through the acquired scores

> in MATLAB

% "scores" is a 2D vector containing the outlier indices provided by ORCA along with the subsequent outlier scores
% "labels" is a vector containing the outlier labels for all data elements; 0 for inliers, and 1 for outliers
% "timElp_Mammography" is the execution time of ORCA on this dataset

scorTmp = zeros(numel(labels),1);
scorTmp(scores(:,1)) = scores(:,2);
scores = scorTmp;

[~,~,~,ROC_Mammography] = perfcurve(labels,scores,1);
[~,~,~,PR_Mammography] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('ORCA result for Mammography:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',ROC_Mammography,PR_Mammography,timElp_Mammography);
save('res_ORCA_Mammography.mat','ROC_Mammography','PR_Mammography','timElp_Mammography');
```


