# EnLOF

## An ensemble version of LOF

EnLOF was firstly introduced in [1] as an ensemble version of LOF and it solely works with the nearest neighbor of each point among the sampled instances in every ensemble member to define the density around it. This method, like _i_NNE, is essentially inspired by the _i_Forest [2][3] method and thus enjoys an adequate number of subsamples (t) with a specific size (![equation](<img src="https://latex.codecogs.com/svg.image?\psi&space;" title="\psi " />)) to determine the anomaly scores for every object. Here, we follow the same premise as _i_Forest and set the two parameters as suggested, i.e., t = 100 and &psi = 256.

[1] Bandaragoda, Tharindu R., et al. "Isolation‐based anomaly detection using nearest‐neighbor ensembles." Computational Intelligence 34.4 (2018): 968-998.

[2] Liu, Fei Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation forest." 2008 eighth ieee international conference on data mining. IEEE, 2008.

[3] Liu, Fei Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation-based anomaly detection." ACM Transactions on Knowledge Discovery from Data (TKDD) 6.1 (2012): 1-39.

## Implementation details

You can follow the subsequent script as a template to use the function and obtain EnLOF results out of an arbitrary dataset:

```matlab
% Setting initial parameters
t = 100; psi = 256; maxIter = 40;

% Mammography dataset
clear X y
load('C:\SDCOR\codes\datasets\realData\Mammography_(11183by6_260o).mat');

tEarr_Mammography = [];
ROCarr_Mammography = [];
PRarr_Mammography = [];
for e1 = 1:maxIter
    tic
    [~,ROC,PR] = EnLOF(X,y,t,psi);
    tEarr_Mammography = [tEarr_Mammography toc];
    ROCarr_Mammography = [ROCarr_Mammography ROC];
    PRarr_Mammography = [PRarr_Mammography PR];
end
ROCavg_Mammography = mean(ROCarr_Mammography); ROCstd_Mammography = std(ROCarr_Mammography);
PRavg_Mammography = mean(PRarr_Mammography); PRstd_Mammography = std(PRarr_Mammography);
timElpAvg_Mammography = mean(tEarr_Mammography);

fprintf('EnLOF (t=%d,psi=%d) result with maxIter = %d for Mammography:\t\tROC = %0.3f+-%0.3f\t\tPR = %0.3f+-%0.3f\t\telpsTime = %0.3f sec\n\n',...
    t,psi,maxIter,ROCavg_Mammography,ROCstd_Mammography,PRavg_Mammography,PRstd_Mammography,timElpAvg_Mammography);
save(['res_EnLOF(t=' num2str(t) ',psi=' num2str(psi) ')_Mammography.mat'],'ROCarr_Mammography','PRarr_Mammography','ROCavg_Mammography','ROCstd_Mammography',...
	'PRavg_Mammography','PRstd_Mammography','timElpAvg_Mammography');
```


