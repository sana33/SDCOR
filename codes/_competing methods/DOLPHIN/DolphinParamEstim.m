function [R] = DolphinParamEstim(X,Eps,delta,alpha,sigma)
% This procedure estimates the suitable value for parameter R
% w.r.t. alpha, the percentage of expected outliers
% 
% X: the query dataset with N objects
% Eps: the estimation error
% delta: the probability constant
% alpha: the expected fraction of outliers to be detected
% sigma: the fraction of neighbors to be considered

% calculating the sample size
n = nComp(alpha,Eps,delta);

% sampling
N = size(X,1);
S = X(randperm(N,n),:);

% locating R
[~,D] = knnsearch(S,S,'K',round(sigma*n)+1);
Dsort = sort(D(:,end),'descend');

cut = floor(alpha*n);
R = mean([Dsort(cut),Dsort(cut+1)]);

end

% Subfunctions here

function [n] = nComp(alpha,Eps,delta)

n = ceil((alpha*(1-alpha)/Eps^2)*norminv(1-delta/2)^2);

end

