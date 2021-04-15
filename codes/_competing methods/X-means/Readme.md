# X-means

## A fast K-means variant optimized for clustering large-scale data

Actually, X-means is not an anomaly recognition technique in essence as it assumes the input data free of noise. Therefore, after obtaining the final clustering outcome through this method, the Euclidean distance of every object to its closest centroid is assigned to it as an outlier score; hence, the following assessment measures could be calculated.

For the evaluations, the minimum and the maximum number of clusters are set to 1 and 15, respectively; the number of times to split a cluster and the maximum number of iterations are set to 6 and 50, respectively, as well. The other interesting parameters for building the KD-tree data structure are *max_leaf_size* and *min_box_width* that are set as suggested by the authors; for the datasets with a cardinality less than 100 thousand points, values equal to 40 and 0.03 are employed, respectively, and for those with greater size, in order, 80 and 0.1 are utilized.

[1] Pelleg, Dan, and Andrew W. Moore. "X-means: Extending k-means with efficient estimation of the number of clusters." Icml. Vol. 1. 2000.

## Implementation details

**The implementation is in C and provided by the genuine authors through [this link](https://www.cs.cmu.edu/~dpelleg/kmeans.html).**





		In the case of datasets in MAT format, you shall convert them to a format acceptable by ORCA. For this matter, first, you should save them through MATLAB in ASCII format; this can be done using the `dlmwrite()` or the `writematrix()` functions. Then you must employ `dprep.exe` given by the ORCA authors along with the mentioned necessities to change the ASCII data into a binary file required by ORCA. Finally, for running the code, it will be only required to use `orca.exe` with the suggested parameters to obtain the anomaly scores. You can utilize the `ptime.exe` executable code to calculate the runtime for each command and save the command output in a specific file by adding `> orca_output.txt` to the end of the command line; outlier scores and the execution time can be elicited out of this file.

You can follow the subsequent script with the suggested parameters as a template to use the ORCA implementation code and obtain the required results out of an arbitrary dataset:

```matlab
%%% Mammography dataset

%% converting the MAT dataset into a binary format acceptable by DOLPHIN

> in MATLAB:

load('Mammography_(11183by6_260o).mat');
DOLPHIN_dssave('mammographyBin',X);
[ds,rows,cols] = DOLPHIN_dsload('mammographyBin'); % just for checking the correctness of the output binary file

%% computing various R values w.r.t. different alpha values

> in MATLAB:

Eps = .01; delta = .1; alpha = [1:10]./100; sigma = .01;
alphaK = numel(alpha);

% Mammography
load('Mammography_(11183by6_260o).mat');
R_Mammography = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_Mammography(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end





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


% mammography
clear
X = load('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\_Competing methods\X-means\datasets\realData\mammography');
y = load('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\_Competing methods\X-means\results\realData\algOutput\mammography_lab');
ctrs = load('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\_Competing methods\X-means\results\realData\algOutput\mammography_ctrs.out');

tic
[~,scores] = knnsearch(ctrs,X);
[~,~,~,ROC_mammography] = perfcurve(y,scores,1);
[~,~,~,PR_mammography] = perfcurve(y,scores,1,'XCrit','reca','YCrit','prec');
timElp_mammography = 1.115+toc;

fprintf('X-means result_02 for mammography:\t\tROC AUC = %0.3f\t\tPR AUC = %0.3f\t\telpsTime = %0.3f sec\n\n',ROC_mammography,PR_mammography,timElp_mammography);
save('res02_X-means_mammography.mat','ROC_mammography','PR_mammography','timElp_mammography');

```


