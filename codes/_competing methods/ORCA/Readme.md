# ORCA

## ORCA (Optimal Reciprocal Collision Avoidance)



## Implementation details

You can follow the subsequent script with the suggested parameters as a template to use the `EnLOF.m` function and obtain the required results out of an arbitrary dataset:

```matlab
% Mammography
scores = load('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\_Competing methods\ORCA\results\realData\Mammography.scores');
labels = load('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\_Competing methods\ORCA\results\realData\Mammography.labels');

tic
scorTmp = zeros(numel(labels),1);
scorTmp(scores(:,1)) = scores(:,2);
scores = scorTmp;

[~,~,~,ROC_Mammography] = perfcurve(labels,scores,1);
[~,~,~,PR_Mammography] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');
timElp_Mammography = 3.353+2.556+toc;

fprintf('ORCA result for Mammography:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',ROC_Mammography,PR_Mammography,timElp_Mammography);
save('res_ORCA_Mammography.mat','ROC_Mammography','PR_Mammography','timElp_Mammography');
```


