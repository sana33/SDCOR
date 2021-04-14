# DOLPHIN

## DOLPHIN (Detecting OutLiers PusHing objects into an INdex)

Angiulli and Fassetti [1] propose DOLPHIN (Detecting OutLiers PusHing objects into an INdex), a distance-based outlier mining method that is particularly dedicated to operating on disk-resident data and functions efficiently in terms of CPU and I/O cost two at a time. Furthermore, both theoretical and practical proofs are presented that the proposed method occupies the memory space as much as a small portion of the dataset. DOLPHIN obtains its efficiency certainly through combining three policies in an integrated plan, namely: 1) careful selection of objects to be retained in the buffer; 2) employing decent pruning strategies; 3) applying effective similarity inspection approaches, which is especially achieved without requiring prior indexing the entire input data, differently from the other competing techniques. Moreover, DOLPHIN is capable of being applied on any kind of records attributed to either metric or non-metric scopes.

For executing the DOLPHIN method on any query data, we need two distinct input parameters specialized for every dataset and two other general parameters that could be established globally in all experiments. The two specific parameters are, namely, $ R $, the neighborhood radius, and $ k $, the minimum neighborhood cardinality required for an object to be identified as an inlier. The other two general parameters are viz \textit{p\textsubscript{inliers}}, the fraction of granted inliers to be maintained in the indexing structure, and \textit{h}, the number of histogram bins used to approximate the nearest neighbors distribution for every point; these parameters are in connection to the pruning rules employed by the method and are determined equal to 0.05 and 16, respectively.

To define $ k $ under each dataset, we have set it to 1\% of the dataset size. However, for $ R $, we followed the \textit{DolphinParamEstim} procedure stipulated in the original paper. Concerning this procedure, the parameter $ R $ directly correlates with the expected ratio of outliers, \textit{alpha}, which is anticipated to be detected by DOLPHIN.

It should be noted that DOLPHIN is the only method in our evaluations that does not provide any anomaly scores for the data elements, and its output is all and solely the definite list of potential outliers. In such a case, the ROC and PR curves will not be appealingly smooth, and the following AUC values will not be very reliable either\footnote{In fact, in more convenient computational conditions, DOLPHIN mostly leads to non-promising detection results. On the other hand, in more compelling parameter settings, i.e., lower values for $ R $ and greater values for $ k $, in spite of higher data-processing costs, the DOLPHIN algorithm outputs incline to be quite deterministic and auspicious in all cases; thus, the subsequent detection accuracy outcomes will be more dependable.}. For this reason, we decided to run the method by various $ R $ values, which are in accordance with different \textit{alpha} values\footnote{In all assessments, \textit{alpha} takes ten different values; 1 to 10 with the step length of 1 for the effectiveness experimentation, and 5 to 50 with the step length of 5 in the case of efficiency test where the outliers ratio in the input data is much higher than usual.}, and rank the detected outliers in the entire iterations w.r.t. the sum of their appearance times in diverse iterations. This heuristic strategy would lead to some sort of outlier ranking, in which every potential outlier gains a positive integer score with a direct relationship to its anomalousness degree, while non-outlier objects attain a score of zero.


[1] Angiulli, Fabrizio, and Fabio Fassetti. "Dolphin: An efficient algorithm for mining distance-based outliers in very large datasets." ACM Transactions on Knowledge Discovery from Data (TKDD) 3.1 (2009): 1-57.

## Implementation details

**The implementation is in C++ and is provided by the genuine authors**. In the case of datasets in MAT format, you shall convert them to a format acceptable by ORCA. For this matter, first, you should save them through MATLAB in ASCII format; this can be done using the `dlmwrite()` or the `writematrix()` functions. Then you must employ `dprep.exe` given by the ORCA authors along with the mentioned necessities to change the ASCII data into a binary file required by ORCA. Finally, for running the code, it will be only required to use `orca.exe` with the suggested parameters to obtain the anomaly scores. You can utilize the `ptime.exe` executable code to calculate the runtime for each command and save the command output in a specific file by adding `> orca_output.txt` to the end of the command line; outlier scores and the execution time can be elicited out of this file.

You can follow the subsequent script with the suggested parameters as a template to use the ORCA implementation code and obtain the required results out of an arbitrary dataset:

```matlab
%%% Mammography dataset

%% converting the MAT dataset into a binary format acceptable by DOLPHIN

(1) in MATLAB:
load('Mammography_(11183by6_260o).mat');
DOLPHIN_dssave('mammographyBin',X);
% [ds,rows,cols] = DOLPHIN_dsload('mammographyBin');



load('Mammography_(11183by6_260o).mat');
dlmwrite('Mammography',X,'precision','%.15f'); dlmwrite('labels',y);

(2) create the "Mammography.fields" file with the following content:
attrib01: continuous.
attrib02: continuous.
attrib03: continuous.
attrib04: continuous.
attrib05: continuous.
attrib06: continuous.

(3) in command window:
dprep.exe Mammography Mammography.fields Mammography.bin -rand -snone -cleanf

%% running the C++ executable code in command window

C:\ptime.exe orca.exe Mammography.bin Mammography.bin Mammography.weights -n 1397 > Mammography_ORCA.comOut

%% gaining the AUC outcomes through the acquired scores

% "scores" contains the outliers indices provided by ORCA along with the subsequent outlier scores
% "labels" contains the outlier labels for all data elements; 0 for inliers, and 1 for outliers
% "timElp_Mammography" is the execution time of ORCA on this dataset

scorTmp = zeros(numel(labels),1);
scorTmp(scores(:,1)) = scores(:,2);
scores = scorTmp;

[~,~,~,ROC_Mammography] = perfcurve(labels,scores,1);
[~,~,~,PR_Mammography] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('ORCA result for Mammography:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',ROC_Mammography,PR_Mammography,timElp_Mammography);
save('res_ORCA_Mammography.mat','ROC_Mammography','PR_Mammography','timElp_Mammography');
```


