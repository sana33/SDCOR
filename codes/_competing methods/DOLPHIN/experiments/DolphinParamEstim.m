function [R] = DolphinParamEstim(X,Eps,delta,alpha,sigma)
% X is the query dataset with N objects
% sigma is the fraction of neighbors to be considered
% alpha is the expected fraction of outliers to be detected
% Eps is the estimation error
% delta is the probability constant

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

