
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ce.aut.ac.ir/~sann_cv/
% June 2020

function [lofVals,lofKmat,AUC_LOF,tElapsed] = LOF(H)

global BLK_SZ_LIM
minPtsLB = H.minPtsIntv(1);
minPtsUB = H.minPtsIntv(end);
kStep = H.kStepLngth;

tStart = tic; % Set start time
matlMat_Maker();

kVals = minPtsLB:kStep:minPtsUB;
kCount = length(kVals);
lofKmat = zeros(n,kCount);

for c1 = 1:kCount
    kDistX = matlMatX(:,kVals(c1));
    kDistNeighbLog = matlMatX <= kDistX;
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

% Set end time
tElapsed = toc(tStart);

% Calculating AUC
[~,~,~,AUC_LOF] = perfcurve(H.labDS.y,lofVals,1);


%% Nested Functions Here!

% Nested function for calculating the neighborhood graph block-by-block
    function matlMat_Maker()
        
        [n,~] = size(H.labDS,'X');
        
        %------- Error handing -------%
        [~,sysView] = memory;
        if n*minPtsUB*8+BLK_SZ_LIM^2*8*6+BLK_SZ_LIM*minPtsUB*8*3 > sysView.PhysicalMemory.Available
            % Creating the structure of the message box
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgCont = {'\fontsize{10}There is not {\bf{enough memory}} to compute the {\bf{materialization matrix}} for LOF!, due to the incorrect input parameter {\bf{minPtsUB}}!'; '\fontsize{10}Thus, please {\color{red}\bf{stop}} the execution and run it again with the true optimal parameters, or even just free up some memory space! ;-)'; '\newline{\color{red}\bf{Note:}} Consider the {\bf{BlckSzLim}} value in GUI as a parameter too!'};
            h = msgbox(msgCont,'Memory Error','error',CreateStruct);
            uiwait(h)
            keyboard % Stop the operation from here, because of lack of memory
        end
        %-----------------------------%

        matlMatX = [];
        matlMatXind = [];
        cntLim = ceil(n/BLK_SZ_LIM);
        
        for d1 = 1:cntLim
            %     d1
            if d1 ~= cntLim
                Y1 = H.labDS.X((d1-1)*BLK_SZ_LIM+1:d1*BLK_SZ_LIM,:);
            else
                Y1 = H.labDS.X((d1-1)*BLK_SZ_LIM+1:end,:);
            end
            chnkSz = size(Y1,1);
            matlMatY = [];
            matlMatYind = [];
            
            for d2 = 1:cntLim
                %         d2
                if d2 ~= cntLim
                    indVec = (d2-1)*BLK_SZ_LIM+1:d2*BLK_SZ_LIM;
                else
                    indVec = (d2-1)*BLK_SZ_LIM+1:n;
                end
                Y2 = H.labDS.X(indVec,:);
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
                
                H.progLevl_statText.String = [num2str((d1-1)*cntLim+d2) '/' num2str(cntLim^2)]; pause(.001);
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

