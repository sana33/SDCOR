function [score,ROC,PR] = EnLOF(X,y,t,psi)

n = size(X,1);
score = ones(n,t);
for c1 = 1:t
    S = X(randperm(n,psi),:);
    
    NSmdl = createns(S);
    
    [~,S1nDst] = knnsearch(NSmdl,S,'K',2);
    S1nDst = S1nDst(:,2);
    [XS1nIdx,XS1nDst] = knnsearch(NSmdl,X);
    
    anmlCond = XS1nDst>S1nDst(XS1nIdx);
    score(anmlCond,c1) = XS1nDst(anmlCond)./S1nDst(XS1nIdx(anmlCond));
end

score = mean(score,2);
[~,~,~,ROC] = perfcurve(y,score,1);
[~,~,~,PR] = perfcurve(y,score,1,'XCrit','reca','YCrit','prec');

end

