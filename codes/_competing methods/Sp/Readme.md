# <img src="https://latex.codecogs.com/svg.image?S_p" title="S_p" />

## Rapid Distance-Based Outlier Detection via Sampling

<img src="https://latex.codecogs.com/svg.image?\inline&space;S_p" title="\inline S_p" /> [1] is a simple and rapid distance-based method that utilizes the nearest neighbor distance on a small sample from the dataset. It takes a small random sample of the entire dataset and then assigns an outlierness score to each point, as the distance from the point to its nearest neighbor in the sample set. Therefore, this method enjoys a linear time complexity concerning each one of the essential variables, viz the number of objects, the number of dimensions, and the number of samples; furthermore, <img src="https://latex.codecogs.com/svg.image?\inline&space;S_p" title="\inline S_p" /> has a constant space complexity, which makes it ideal for analyzing massive datasets.

<img src="https://latex.codecogs.com/svg.image?\inline&space;S_p" title="\inline S_p" /> is really simple to understand and also to implement, and only requires one single parameter, which even with its default value proposed by the authors, promising outcomes could be achieved over multiple datasets. For this reason, we follow the same procedure as suggested in the original paper and set the sample size, _s_, equal to 20 in our experiments.

[1] Sugiyama, Mahito, and Karsten Borgwardt. "Rapid distance-based outlier detection via sampling." Advances in Neural Information Processing Systems 26 (2013): 467-475.

## Implementation details

You can follow the subsequent script with the suggested parameters as a template to use the `Sp.m` function and obtain the required results out of an arbitrary dataset:

```matlab
% Setting initial parameters
s = 20; % sample size
totIter = 40; % total number of independent runs

% Mammography dataset
load('Mammography_(11183by6_260o).mat');

tEarr_Mammography = [];
ROCarr_Mammography = [];
PRarr_Mammography = [];
for e1 = 1:totIter
    tic
    [~,ROC,PR] = Sp(X,y,s);
    tEarr_Mammography = [tEarr_Mammography toc];
    ROCarr_Mammography = [ROCarr_Mammography ROC];
    PRarr_Mammography = [PRarr_Mammography PR];
end
ROCavg_Mammography = mean(ROCarr_Mammography); ROCstd_Mammography = std(ROCarr_Mammography);
PRavg_Mammography = mean(PRarr_Mammography); PRstd_Mammography = std(PRarr_Mammography);
timElpAvg_Mammography = mean(tEarr_Mammography);

fprintf('Sp (s=%d) result with totIter = %d for Mammography:\t\tROC = %0.3f+-%0.3f\t\tPR = %0.3f+-%0.3f\t\telpsTime = %0.3f sec\n\n',...
    s,totIter,ROCavg_Mammography,ROCstd_Mammography,PRavg_Mammography,PRstd_Mammography,timElpAvg_Mammography);
save(['res_Sp(s=' num2str(s) ')_Mammography.mat'],'ROCarr_Mammography','PRarr_Mammography','ROCavg_Mammography','ROCstd_Mammography','PRavg_Mammography','PRstd_Mammography','timElpAvg_Mammography');
```
