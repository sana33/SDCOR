function [R,ROC,PR] = Sp(X,y,s)
% % This is an implementation for the distance-based method for outlier
% % detection named Sp, regarding the paper "Rapid Distance-Based Outlier
% % Detection via Sampling", by Mahito Sugiyama and Karsten M. Borgwardt
% 
% X: input data set
% y: label vector
% s: sample size
 

n = size(X,1);
S = X(randperm(n,s),:);
[~,R] = knnsearch(S,X);

[~,~,~,ROC] = perfcurve(y,R,1);
[~,~,~,PR] = perfcurve(y,R,1,'XCrit','reca','YCrit','prec');

end

