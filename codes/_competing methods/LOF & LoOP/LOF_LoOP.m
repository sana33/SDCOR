function [lofVals,lofKmat,ROC_LOF,PR_LOF,tElapsed_LOF,LoOPvals,ROC_LoOP,PR_LoOP,tElapsed_LoOP] = LOF_LoOP(fileName,X,y,minPtsLB,minPtsUB,kStep,...
    kVal,lambda,OlNbCnd,NSMethod,blkSzLim)

% minPtsLB, minPtsUB, and kStep are LOF parameters
% kVal and lambda are LoOP parameters
% OlNbCnd is the Overlapping Neighborhood (OlNb) Condition for the materialization matrix
% NSMethod is the Neighbor Search Method for the knnsearch function [0:kdtree, 1:exhaustive] while OlNb is equal to 0
% blkSzLim is the memory size limit while OlNb equals 1 (our very own implementation)

if ~OlNbCnd
    switch NSMethod
        case 0
            OlNbMod = 'knn-nonOlNb_NSMeth=kdtree';
        case 1
            OlNbMod = 'knn-nonOlNb_NSMeth=exhaustive';
    end
else
    OlNbMod = 'OlNb';
end
n = size(X,1);

tic % Set start time for the mateialization matrix
matlMat_Maker();
matlMat_time = toc; % Set end time for the mateialization matrix

tic % Set start time for LOF
kVals = minPtsLB:kStep:minPtsUB;
kCount = length(kVals);
lofKmat = zeros(n,kCount);

for c1 = 1:kCount
    if OlNbCnd
        kDistX = matlMatX(:,kVals(c1));
        kDistNeighbLog = matlMatX <= kDistX;
    else
        kDistNeighbLog = [ones(n,kVals(c1)) zeros(n,minPtsUB-kVals(c1))];
    end
    kDistNghbCard = sum(kDistNeighbLog,2);
    reachDistMat = max(reshape(matlMatX(matlMatXind(:),kVals(c1)),size(matlMatX)),matlMatX);
    reachDistMatSum = reachDistMat.*kDistNeighbLog;
    reachDistMatSum(isnan(reachDistMatSum) | isinf(reachDistMatSum)) = 0;
    lrdK = kDistNghbCard./sum(reachDistMatSum,2);
    lrdK(isnan(lrdK) | isinf(lrdK)) = 0;
    lrdKmat = reshape(lrdK(matlMatXind(:)),size(matlMatX));
    lofKmat(:,c1) = sum(lrdKmat.*kDistNeighbLog,2)./lrdK./kDistNghbCard;
end

lofKmat(isnan(lofKmat) | isinf(lofKmat)) = 0;

lofVals = max(lofKmat,[],2);

tElapsed_LOF = toc+matlMat_time; % Set end time for LOF

% Calculating ROC and PR for LOF
[~,~,~,ROC_LOF] = perfcurve(y,lofVals,1);
[~,~,~,PR_LOF] = perfcurve(y,lofVals,1,'XCrit','reca','YCrit','prec');

% Printing and saving the output of LOF
fprintf('\nROC = %0.3f and PR = %0.3f obtained out of LOF_%s (MinPtsIntv = [%d,%d], kStep = %d) in %0.3f secs for dataset %s\n',...
    ROC_LOF,PR_LOF,OlNbMod,minPtsLB,minPtsUB,kStep,tElapsed_LOF,fileName);
save(sprintf('res_LOF_%s_(MinPtsIntv=[%d,%d],kStep=%d)_%s_ROC=%0.3f_PR=%0.3f_time=%0.3f.mat',OlNbMod,minPtsLB,minPtsUB,kStep,...
    fileName,ROC_LOF,PR_LOF,tElapsed_LOF),'lofVals','ROC_LOF','PR_LOF','tElapsed_LOF','minPtsLB','minPtsUB','kStep');

%% LoOP operations
tic % Set start time for LoOP
if OlNbCnd
    kDistX = matlMatX(:,kVal);
    kDistNeighbLog = matlMatX <= kDistX;
else
    kDistNeighbLog = [ones(n,kVals(c1)) zeros(n,minPtsUB-kVals(c1))];
end
kDistNghbCard = sum(kDistNeighbLog,2);

Sdist = matlMatX.*kDistNeighbLog; Sdist(isnan(Sdist) | isinf(Sdist)) = 0;

sigma = sqrt(sum(Sdist.^2,2)./kVal);

pDist = lambda*sigma;

pDistNghb = pDist(matlMatXind).*kDistNeighbLog; pDistNghb(isnan(pDistNghb) | isinf(pDistNghb)) = 0;
Epdist = sum(pDistNghb,2)./kDistNghbCard;

pdEpd = pDist./Epdist; pdEpd(isnan(pdEpd) | isinf(pdEpd)) = 0;
plof = pdEpd - 1;

nPLOF = lambda*sqrt(mean(plof.^2));

erfContent = plof./(nPLOF * sqrt(2)); erfContent(isnan(erfContent) | isinf(erfContent)) = 0;
LoOPvals = erf(erfContent);
LoOPvals(LoOPvals<0) = 0;

tElapsed_LoOP = toc+matlMat_time; % Set end time for LoOP

% Calculating ROC and PR for LoOP
[~,~,~,ROC_LoOP] = perfcurve(y,LoOPvals,1);
[~,~,~,PR_LoOP] = perfcurve(y,LoOPvals,1,'XCrit','reca','YCrit','prec');

% Printing the output of the algorithm
fprintf('\nROC = %0.3f and PR = %0.3f obtained out of LoOP_%s (kVal = %d, lambda = %d) in %0.3f secs for dataset %s\n\n',...
    ROC_LoOP,PR_LoOP,OlNbMod,kVal,lambda,tElapsed_LoOP,fileName);
save(sprintf('res_LoOP_%s_(kVal=%d,lambda=%d)_%s_ROC=%0.3f_PR=%0.3f_time=%0.3f.mat',OlNbMod,kVal,lambda,fileName,ROC_LoOP,...
    PR_LoOP,tElapsed_LoOP),'LoOPvals','ROC_LoOP','PR_LoOP','tElapsed_LoOP','kVal','lambda');


%% Nested Functions Here!

% Nested function for calculating the neighborhood graph block-by-block
    function matlMat_Maker()
        
        if ~OlNbCnd
            %------- Error handing -------%
            [~,sysView] = memory;
            if n*minPtsUB*8 > sysView.PhysicalMemory.Available
                % Creating the structure of the message box
                CreateStruct.Interpreter = 'tex';
                CreateStruct.WindowStyle = 'modal';
                msgCont = {'\fontsize{10}There is not {\bf{enough memory}} to compute the {\bf{materialization matrix}} for LOF/LoOP!, due to the incorrect input parameter {\bf{minPtsUB}}!'; '\fontsize{10}Thus, please {\color{red}\bf{stop}} the execution and run it again with the true optimal parameters, or even just free up some memory space! ;-)';};
                h = msgbox(msgCont,'Memory Error','error',CreateStruct);
                uiwait(h)
                keyboard % Stop the operation from here, because of the lack of memory
            end
            %-----------------------------%
            
            switch NSMethod
                case 0
                    [matlMatXind,matlMatX] = knnsearch(X,X,'K',minPtsUB+1,'NSMethod','kdtree');
					matlMatXind = matlMatXind(:,2:end); matlMatX = matlMatX(:,2:end);
                case 1
                    [matlMatXind,matlMatX] = knnsearch(X,X,'K',minPtsUB+1,'NSMethod','exhaustive');
					matlMatXind = matlMatXind(:,2:end); matlMatX = matlMatX(:,2:end);
            end
            return
        end
        
        %------- Error handing -------%
        [~,sysView] = memory;
        if n*minPtsUB*8+blkSzLim^2*8*6+blkSzLim*minPtsUB*8*3 > sysView.PhysicalMemory.Available
            % Creating the structure of the message box
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgCont = {'\fontsize{10}There is not {\bf{enough memory}} to compute the {\bf{materialization matrix}} for LOF/LoOP!, due to the incorrect input parameter {\bf{minPtsUB}}!'; '\fontsize{10}Thus, please {\color{red}\bf{stop}} the execution and run it again with the true optimal parameters, or even just free up some memory space! ;-)'; '\newline{\color{red}\bf{Note:}} Consider the {\bf{BlckSzLim}} value in GUI as a parameter too!'};
            h = msgbox(msgCont,'Memory Error','error',CreateStruct);
            uiwait(h)
            keyboard % Stop the operation from here, because of the lack of memory
        end
        %-----------------------------%
        
        matlMatX = [];
        matlMatXind = [];
        cntLim = ceil(n/blkSzLim);
        
        for d1 = 1:cntLim
            %     d1
            if d1 ~= cntLim
                Y1 = X((d1-1)*blkSzLim+1:d1*blkSzLim,:);
            else
                Y1 = X((d1-1)*blkSzLim+1:end,:);
            end
            chnkSz = size(Y1,1);
            matlMatY = [];
            matlMatYind = [];
            
            for d2 = 1:cntLim
                %         d2
                if d2 ~= cntLim
                    indVec = (d2-1)*blkSzLim+1:d2*blkSzLim;
                else
                    indVec = (d2-1)*blkSzLim+1:n;
                end
                Y2 = X(indVec,:);
                [sDistY, sDistYind] = sort([pdist2(Y1,Y2) matlMatY],2);
                
                if d1 ~= d2
                    kDist = sDistY(:,minPtsUB);
                    kNNlastInd = find(any(sDistY<=kDist),1,'last');
                    matlMatY = sDistY(:,1:kNNlastInd);
                    
                    truInd = [repmat(indVec,chnkSz,1) matlMatYind];
                    sDistYlinInd = (sDistYind-1).*chnkSz+[1:chnkSz]';
                    matlMatYind = truInd(sDistYlinInd);
                    matlMatYind = matlMatYind(:,1:kNNlastInd);
                else
                    kDist = sDistY(:,minPtsUB+1);
                    kNNlastInd = find(any(sDistY<=kDist),1,'last');
                    matlMatY = sDistY(:,2:kNNlastInd);
                    
                    truInd = [repmat(indVec,chnkSz,1) matlMatYind];
                    sDistYlinInd = (sDistYind-1).*chnkSz+[1:chnkSz]';
                    matlMatYind = truInd(sDistYlinInd);
                    matlMatYind = matlMatYind(:,2:kNNlastInd);
                end
                
                % Printing progress of building the materialization matrix
                fprintf('iteration %d of total %d, for dataset %s\n',cntLim*(d1-1)+d2,cntLim^2,fileName);
            end
            
            Xwdth = size(matlMatX,2);
            Ywdth = size(matlMatY,2);
            if Xwdth<Ywdth
                matlMatX = [matlMatX inf*ones(size(matlMatX,1),Ywdth-Xwdth); matlMatY];
                matlMatXind = [matlMatXind n*ones(size(matlMatXind,1),Ywdth-Xwdth); matlMatYind];
            elseif Xwdth>Ywdth
                matlMatX = [matlMatX; matlMatY inf*ones(size(matlMatY,1),Xwdth-Ywdth)];
                matlMatXind = [matlMatXind; matlMatYind n*ones(size(matlMatYind,1),Xwdth-Ywdth)];
            else
                matlMatX = [matlMatX; matlMatY];
                matlMatXind = [matlMatXind; matlMatYind];
            end
        end
        
    end

end

