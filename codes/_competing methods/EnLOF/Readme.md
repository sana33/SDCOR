# EnLOF

## An ensemble version of LOF

EnLOF was firstly introduced in [1] as an ensemble version of LOF and it solely works with each point nearest neighbor among the sampled instances in every ensemble member to define the density around it and the resultant anomaly score. This method, like _i_-NNE [1], is essentially inspired by the _i_-Forest method [2][3] and thus enjoys an adequate number of subsamples (<img src="https://latex.codecogs.com/svg.image?t" title="t" />) with a specific size (<img src="https://latex.codecogs.com/svg.image?\psi&space;" title="\psi " />) to determine the anomaly scores for every object. Here, we follow the same premise as _i_-Forest and set the two parameters as suggested, i.e., <img src="https://latex.codecogs.com/svg.image?t&space;=&space;100" title="t = 100" /> and <img src="https://latex.codecogs.com/svg.image?\psi&space;=&space;256" title="\psi = 256" />.

[1] Bandaragoda, Tharindu R., et al. "Isolation-based anomaly detection using nearest-neighbor ensembles." Computational Intelligence 34.4 (2018): 968-998.

[2] Liu, Fei Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation forest." 2008 eighth ieee international conference on data mining. IEEE, 2008.

[3] Liu, Fei Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation-based anomaly detection." ACM Transactions on Knowledge Discovery from Data (TKDD) 6.1 (2012): 1-39.

## Implementation details

You can follow the subsequent script with the suggested parameters as a template to use the function and obtain EnLOF results out of an arbitrary dataset:

```matlab
% Setting initial parameters
t = 100; % ensemble size
psi = 256; % subsample size
totIter = 40; % total number of independent runs

% Mammography dataset
clear X y
load('C:\SDCOR\codes\datasets\realData\Mammography_(11183by6_260o).mat');

tEarr_Mammography = [];
ROCarr_Mammography = [];
PRarr_Mammography = [];
for e1 = 1:totIter
    tic
    [~,ROC,PR] = EnLOF(X,y,t,psi);
    tEarr_Mammography = [tEarr_Mammography toc];
    ROCarr_Mammography = [ROCarr_Mammography ROC];
    PRarr_Mammography = [PRarr_Mammography PR];
end
ROCavg_Mammography = mean(ROCarr_Mammography); ROCstd_Mammography = std(ROCarr_Mammography);
PRavg_Mammography = mean(PRarr_Mammography); PRstd_Mammography = std(PRarr_Mammography);
timElpAvg_Mammography = mean(tEarr_Mammography);

fprintf('EnLOF (t=%d,psi=%d) result with totIter = %d for Mammography:\t\tROC = %0.3f+-%0.3f\t\tPR = %0.3f+-%0.3f\t\telpsTime = %0.3f sec\n\n',...
    t,psi,totIter,ROCavg_Mammography,ROCstd_Mammography,PRavg_Mammography,PRstd_Mammography,timElpAvg_Mammography);
save(['res_EnLOF(t=' num2str(t) ',psi=' num2str(psi) ')_Mammography.mat'],'ROCarr_Mammography','PRarr_Mammography','ROCavg_Mammography','ROCstd_Mammography',...
	'PRavg_Mammography','PRstd_Mammography','timElpAvg_Mammography');
```


