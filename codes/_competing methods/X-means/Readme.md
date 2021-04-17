# X-means

## A fast K-means variant optimized for clustering large-scale data

Actually, X-means is not an anomaly recognition technique in essence as it assumes the input data free of noise. Therefore, after obtaining the final clustering outcome through this method, the Euclidean distance of every object to its closest centroid is assigned to it as an outlier score; hence, the following assessment measures could be calculated.

For the evaluations, the minimum and the maximum number of clusters are set to 1 and 15, respectively; the number of times to split a cluster and the maximum number of iterations are set to 6 and 50, respectively, as well. The other interesting parameters for building the KD-tree data structure are *max_leaf_size* and *min_box_width* that are set as suggested by the authors; for the datasets with a cardinality less than 100 thousand points, values equal to 40 and 0.03 are employed, respectively, and for those with greater size, in order, 80 and 0.1 are utilized.

[1] Pelleg, Dan, and Andrew W. Moore. "X-means: Extending k-means with efficient estimation of the number of clusters." Icml. Vol. 1. 2000.

## Implementation details

**The implementation is in C and provided by the genuine authors through [this link](https://www.cs.cmu.edu/~dpelleg/kmeans.html).**

In the case of datasets in MAT format, you shall convert them to a format acceptable by X-means. For this matter, you should save them through MATLAB in ASCII format; this can be done using the `dlmwrite()` or the `writematrix()` functions. Moreover, the corresponding outlier labels could be saved separately in a file with ASCII format for further accuracy assessments.

Finally, for running the code, it will only required to use the `kmeans.exe` function in the installation directory of X-means. You can find the necessary information on how to employ this function with the suggested parameters in the `README` file in the same directory. Besides, you can utilize the `ptime.exe` executable code to calculate the runtime of every command.

The output of the X-means method is a file containing the final means located by the method. You should utilize this file for calculating the final detection results.

You can follow the subsequent scripts with the suggested parameters as a template to use the X-means implementation code and obtain the required results out of an arbitrary dataset:

```matlab
%%% Mammography dataset
```

```matlab
%% converting the MAT dataset into ASCII format acceptable by X-means

> in MATLAB:

load('Mammography_(11183by6_260o).mat');
dlmwrite('Mammography',X,'precision','%.15f'); dlmwrite('Mammography_lab',y);
```

```matlab
%% running the C executable code of X-means in the command window

> in the command window:

ptime.exe kmeans.exe kmeans -k 1 -create_universe true -D_SHOW_END_CENTERS -method blacklist -max_leaf_size 40 -min_box_width 0.03 -cutoff_factor 0.5 -max_iter 50 -num_splits 6 -max_ctrs 15 -in Mammography -save_ctrs Mammography_ctrs
```

```matlab
%% gaining the AUC outcomes through the acquired cluster centers

> in MATLAB:

% "Mammography" represents the Mammography dataset
% "Mammography_lab" contains the outlier labels for all data elements; 0 for inliers, and 1 for outliers
% "Mammography_ctrs" contains the final cluster centers provided by X-means as output
% "timElp_Mammography" is the execution time of X-means on this dataset

[~,scores] = knnsearch(Mammography_ctrs,Mammography);
[~,~,~,ROC_Mammography] = perfcurve(y,scores,1);
[~,~,~,PR_Mammography] = perfcurve(y,scores,1,'XCrit','reca','YCrit','prec');

fprintf('X-means result for Mammography:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',ROC_Mammography,PR_Mammography,timElp_Mammography);
save('res_X-means_Mammography.mat','ROC_Mammography','PR_Mammography','timElp_Mammography');
```
