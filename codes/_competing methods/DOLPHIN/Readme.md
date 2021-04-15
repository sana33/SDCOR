# DOLPHIN

## DOLPHIN (Detecting OutLiers PusHing objects into an INdex)

Angiulli and Fassetti [1] propose DOLPHIN (Detecting OutLiers PusHing objects into an INdex), a distance-based outlier mining method that is particularly dedicated to operating on disk-resident data and functions efficiently in terms of CPU and I/O cost two at a time. Furthermore, both theoretical and practical proofs are presented that the proposed method occupies the memory space as much as a small portion of the dataset. DOLPHIN obtains its efficiency certainly through combining three policies in an integrated plan, namely: 1) careful selection of objects to be retained in the buffer; 2) employing decent pruning strategies; 3) applying effective similarity inspection approaches, which is especially achieved without requiring prior indexing the entire input data, differently from the other competing techniques. Moreover, DOLPHIN is capable of being applied on any kind of records attributed to either metric or non-metric scopes.

For executing the DOLPHIN method on any query data, we need two distinct input parameters specialized for every dataset and two other general parameters that could be established globally in all experiments. The two specific parameters are, namely, _R_, the neighborhood radius, and _k_, the minimum neighborhood cardinality required for an object to be identified as an inlier. The other two general parameters are viz <img src="https://latex.codecogs.com/svg.image?p_{inliers}" title="p_{inliers}" />, the fraction of granted inliers to be maintained in the indexing structure, and _h_, the number of histogram bins used to approximate the nearest neighbors distribution for every point; these parameters are in connection to the pruning rules employed by the method and are determined equal to 0.05 and 16, respectively.

To define _k_ under each dataset, we have set it to 1% of the dataset size. However, for _R_, we followed the _DolphinParamEstim_ procedure stipulated in the original paper. Concerning this procedure, the parameter _R_ directly correlates with the expected ratio of outliers, _alpha_, which is anticipated to be detected by DOLPHIN.

It should be noted that DOLPHIN does not provide any anomaly scores for the data elements, and its output is all and solely the definite list of potential outliers. In such a case, the ROC and PR curves will not be appealingly smooth, and the following AUC values will not be very reliable either <sup>1</sup>. For this reason, as proposed by the authors, we decided to run the method by various _R_ values, which are in accordance with different _alpha_ values, and rank the detected outliers in the entire iterations w.r.t. the sum of their appearance times in diverse iterations. This heuristic strategy would lead to some sort of outlier ranking, in which every potential outlier gains a positive integer score with a direct relationship to its anomalousness degree, while non-outlier objects attain a score of zero.

<sup>1</sup> *In fact, in more convenient computational conditions, DOLPHIN mostly leads to non-promising detection results. On the other hand, in more compelling parameter settings, i.e., lower values for _R_ and greater values for _k_, in spite of higher data-processing costs, the DOLPHIN algorithm outputs incline to be quite deterministic and auspicious in all cases; thus, the subsequent detection accuracy outcomes will be more dependable.*

[1] Angiulli, Fabrizio, and Fabio Fassetti. "Dolphin: An efficient algorithm for mining distance-based outliers in very large datasets." ACM Transactions on Knowledge Discovery from Data (TKDD) 3.1 (2009): 1-57.

## Implementation details

**The implementation code is a Linux C++ executable binary file provided by the genuine authors.**

In the case of datasets in MAT format, you shall convert them to a format acceptable by DOLPHIN. This can be done through MATLAB codes dedicated to this issue by the authors.

After preparing the input dataset, it is required to compute various _R_ quantities regarding diverse _alpha_ values. This could be achieved through the _DolphinParamEstim_ procedure with the suggested parameters.

Then for running the DOLPHIN method on a query dataset, you should do as the following template in the Linux Terminal window:

	./dolphin <File.ds> <k> <R> <silent (t/f)> <prob> <slots>

	where:

	- <File.ds> is the dataset binary file

	- <k> is the numer of nearest neigbhors to consider

	- <R> is the radius value to consider

	- <silent (t/f)> is a flag to disable(t)/enable(f) verbose mode

	- <prob> is the p_inliers parameter (the suggested value is equal to 0.05)

	- <slots> is the number of histogram bins used to approximate nearest neighbors distribution (the suggested value is equal to 16)

This command will produce some files as output which the one starting with "Outliers_" contains the outlier indices and should be utilized to compute detection accuracy results.

After obtaining outlier indices w.r.t. various _R_ values, we need to acquire AUC results out of them. If there is only one run, you can utilize the `Experiment_singAUCcalc.m` procedure to calculate the AUROC and AUPRC results; in the case of multiple runs, the `Experiment_multAUCcalc.m` procedure could be of use. Note that for outlier labels, you need to extract them from the corresponding dataset file and save it, e.g., in an ASCII format with an extension like LABELS. Moreover, in both procedures, we load the command output file(s) to extract the related execution time(s).

You can follow the subsequent scripts with the suggested parameters as a template to use the DOLPHIN implementation code and obtain the required results out of an arbitrary dataset:

```matlab
%%% Mammography dataset
```

```matlab
%% converting the MAT dataset into a binary format acceptable by DOLPHIN

% Note: As for the dataset format, it must be a binary file of float32 numbers containing data points in row major order (n*d*4 bytes, n=number of rows, d=number of columns). At the beginning of the file, a header is required (8 bytes), consisting of the number of columns (d) and rows (n), respectively, stored as two int32 numbers.

> in MATLAB:

load('Mammography_(11183by6_260o).mat');
DOLPHIN_dssave('mammographyBin',X);
[ds,rows,cols] = DOLPHIN_dsload('mammographyBin'); % just for checking the correctness of the output binary file
```

```matlab
%% computing various R values w.r.t. different alpha values

> in MATLAB:

Eps = .01; delta = .1; alpha = [1:10]./100; sigma = .01;
alphaK = numel(alpha);

load('Mammography_(11183by6_260o).mat');
R_Mammography = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_Mammography(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end
```

```matlab
%% running the Linux C++ binary code of DOLPHIN with k=112 and R=1.170412063

> in Linux Terminal window:

time sudo ./dolphin.bin mammographyBin 112 1.170412063 t 0.05 16 |& tee output_mammography_k=112_r=1.170412063.txt

% Note: "time" and "sudo" commands are for calculating the runtime and giving administrative privileges to the command. "|& tee" is for controlling the command output. "output_mammography_k=112_r=1.170412063.txt" file will contain the command output, including the runtime.
```
