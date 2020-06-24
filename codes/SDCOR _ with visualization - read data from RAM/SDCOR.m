
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ceit.aut.ac.ir/~sann_cv/
% June 2020

function [] = SDCOR(hO,H)

SDCORstrt = tic; % Execution time for SDCOR

global DBDS CLUSTS CLUST_MODF_ARR CHNK_MEMB_COND

timWast = 0; % Calculating the wasted time for visualization
% tSt = tic; % Setting start time for sampling phase @-- debugging script --@
while true
    rndSmp_Maker();
    
    sampPrmChosMthd();
    
    [H.idxSamp,singChck] = posSemiDefCheck(H.sampData,H.idxSamp);
    
    if ~singChck
        break;
    else
        fprintf('Singularity happened during "Sampling" phase! The procedure will be conducted again\n');
        beep
    end
end
% tElp0 = toc(tSt); % Setting elapsed time for sampling phase @-- debugging script --@
% fprintf('\n\nRandSamp:\t%0.2fe-4 sec\n\n',tElp0*1e4); % @-- debugging script --@

tic
if H.dispOn
    plotOptional(H,{},'sampDS');
end
timWast = timWast+toc; % Adding up the wasted time for visualization

CLUSTS = cell(0);
clustInfMaker(H,{H.sampData,H.idxSamp});
H.origK = size(CLUSTS,2);
H.origK_val_statText.String = num2str(H.origK);

H.sampDetArr = cell2mat(CLUSTS(8,:));
H.sampDetArr(H.sampDetArr<=0) = 1; % Error handling

retainSet = [];
retIdx = [];
accResArr = [];
H.accFinCond = 0;

% % Making the array of the number and ratio of sustained objects in RAM belonging to each chunk @-- debugging script --@
% sustObjOfCurrChnkArr = []; % @-- debugging script --@

maxIter = ceil(H.n/H.chunkSz);
for c1 = 1:maxIter
    indStart = (c1-1)*H.chunkSz+1;
    if c1~=maxIter
        indEnd = c1*H.chunkSz;
    else
        indEnd = H.n;
    end
    chunk = H.DS(indStart:indEnd,:);
    
%     tSt = tic; % Setting start time for chunk memb. @-- debugging script --@
    CHNK_MEMB_COND = 1;
    CLUST_MODF_ARR = transpose(1:size(CLUSTS,2));
    [clustBestArr] = clustAccMahal(H,chunk);
    CHNK_MEMB_COND = 0;
    
    clustInfUpdate(chunk,clustBestArr,H.PCvarRat);
%     tElp1 = toc(tSt); % Setting elapsed time for chunk memb. @-- debugging script --@
    
    retainSet = [retainSet; chunk(clustBestArr==0,:)];
    retIdx = [retIdx; (c1-1)*H.chunkSz+find(clustBestArr==0)];
    
%     tElp2 = 0; tElp3 = 0; tElp4 = 0; % @-- debugging script --@
    if ~isempty(retIdx)
%         tSt = tic; % Setting start time for retain set memb. 1 @-- debugging script --@
        [retainSet,retIdx] = retSetClustMembCheck(H,retainSet,retIdx);
%         tElp2 = toc(tSt); % Setting elalsed time for retain set memb. 1 @-- debugging script --@
        
		if ~isempty(retIdx)
% 	        tSt = tic; % Setting start time for retain set clust. @-- debugging script --@
			befrLastCretClstNoK = size(CLUSTS,2);
			DBDS = retainSet;
			[retainSet,retIdx] = retSetClustMaker(H,retainSet,retIdx);
            CLUST_MODF_ARR = transpose(befrLastCretClstNoK+1:size(CLUSTS,2));
% 	        tElp3 = toc(tSt); % Setting elapsed time for retain set clust. @-- debugging script --@
			
            if ~isempty(retIdx)
%                 tSt = tic; % Setting start time for retain set memb. 2 @-- debugging script --@
                [retainSet,retIdx] = retSetClustMembCheck(H,retainSet,retIdx);
%                 tElp4 = toc(tSt); % Setting elalsed time for retain set memb. 2 @-- debugging script --@
            end
		end
        
    end
	
%     fprintf('ChunkMemb:\t%0.2fe-4 sec;\tretSetMemb.1:\t%0.2fe-4 sec;\tretSetClust.:\t%0.2fe-4 sec;\tretSetMemb.2:\t%0.2fe-4 sec\n',...
%         tElp1*1e4,tElp2*1e4,tElp3*1e4,tElp4*1e4); % Printing elapsed times for various stages @-- debugging script --@
    
%     % Demonstrating the number and ratio of sustained objects in RAM belonging to current chunk @-- debugging script --@
%     sustObjNoRAM(1); % @-- debugging script --@
    
    tic
    [accResArr] = accReport(H,accResArr,retIdx);
    timWast = timWast+toc; % Adding up the wasted time for visualization
end

% % Demonstrating the mean value of numbers and ratios of sustained objects in RAM belonging to all processed chunks @-- debugging script --@
% sustObjNoRAM(2); % @-- debugging script --@

tic
H.accFinCond = 1;
H.accResArr = accResArr;
plotOptional(H,{accResArr},'accPerChunk');
timWast = timWast+toc; % Adding up the wasted time for visualization

H.clusts = CLUSTS;
H.retIdx = retIdx;

H.means = cell2mat(transpose(CLUSTS(1,:)));

tic
retSetDisp();
timWast = timWast+toc; % Adding up the wasted time for visualization

% tSt = tic; % Setting start time for building the final clustering model @-- debugging script --@
[H.idxMeans,H.finalClusts,H.meansMeans,H.regenDS,H.idxRegenDS] = finalClustsMaker(H);
% tElp5 = toc(tSt); % Setting elapsed time for building the final clustering model @-- debugging script --@
% fprintf('FinClstMakr:\t%0.2fe-4 sec\n',tElp5*1e4); % @-- debugging script --@

tic
finMnsRegDS_Disp();
timWast = timWast+toc; % Adding up the wasted time for visualization

% tSt = tic; % Setting start time for the scoring phase @-- debugging script --@
[H.mahalScores,H.idxFin,H.finalAUC] = OLscoreMaker(H);
% tElp6 = toc(tSt); % Setting elapsed time for the scoring phase @-- debugging script --@
% fprintf('OLscrMakr:\t%0.2fe-4 sec\n',tElp6*1e4); % @-- debugging script --@

tic
if H.dispOn
    plotOptional(H,{},'scorDS');
end
timWast = timWast+toc; % Adding up the wasted time for visualization


%% Nested functions here!

% Nested function for random sampling
    function rndSmp_Maker()
        H.sampIdx = transpose(randperm(H.n,floor(H.sampRate*H.n)));
        H.sampData = H.DS(H.sampIdx,:);
        
        if H.dispOn
            if H.p>2
                [~,sampData_PCA,~] = pca(H.sampData);
                H.sampData_PCA = sampData_PCA(:,1:2);
            else
                H.sampData_PCA = [];
            end
        else
            H.sampData_PCA = [];
        end
        
    end

% Nested function for selecting the type of DBSCAN parameter choosing
    function sampPrmChosMthd()
        switch H.PCM
            case 'PSO_pcm_radioBtn'
                DBDS = H.sampData;
                [bestSolParams,bestSolIdx,PSO_costArr] = PSO_DBSCAN(H);
                
                H.paramSampDS = bestSolParams;
                H.idxSamp = bestSolIdx;
                H.paramCostArrSamp = {PSO_costArr};
                
                H.epsilonFin = H.epsCoeff*H.paramSampDS(1);
                H.MinPtsFin = ceil(H.MinPtsCoeff*H.paramSampDS(2));
                
            case 'manu_pcm_radioBtn'
                DBDS = H.sampData;
                [H.idxSamp,~] = DBSCAN(H.manuEps,H.manuMnPt);
                
                H.paramSampDS = [H.manuEps H.manuMnPt];
                H.paramCostArrSamp = {};
                
                H.epsilonFin = H.epsCoeff*H.manuEps;
                H.MinPtsFin = ceil(H.MinPtsCoeff*H.manuMnPt(1));
        end
    end

% Nested function for displaying retain set
    function retSetDisp()
        if H.dispOn
            if H.p>2
                H.means_PCA = H.means*H.coef_PCA;
            else
                H.means_PCA = [];
            end
            plotOptional(H,{},'retSetDS');
        else
            H.means_PCA = [];
        end
    end

% Nested function for displayign final means
    function finMnsRegDS_Disp()
        if H.dispOn
            if H.p>2
                H.meansMeans_PCA = H.meansMeans*H.coef_PCA;
                regenDS_PCA = H.regenDS*H.coef_PCA;
                H.regenDS_PCA = regenDS_PCA(:,1:2);
            else
                H.meansMeans_PCA = [];
                H.regenDS_PCA = [];
            end
            H.origKvec = transpose(1:H.origK);
            plotOptional(H,{},'finalMeans');
        else
            H.origKvec = [];
            H.meansMeans_PCA = [];
            H.regenDS_PCA = [];
        end
        
        if H.dispOn
            plotOptional(H,{},'regenDS');
        end
    end

% Nested function for displaying the number of sustained objects of current chunk
    function sustObjNoRAM(iterType)
        switch iterType
            case 1
                sustObjOfCurrChnk = sum(clustBestArr==0);
                sustObjOfCurrChnkArr = [sustObjOfCurrChnkArr [sustObjOfCurrChnk; sustObjOfCurrChnk/H.chunkSz*100]];
                fprintf('No. and ratio of sustained objects of chunk %d is:\t%d\t%0.2f%%\n',c1,sustObjOfCurrChnkArr(1,end), ...
                    sustObjOfCurrChnkArr(2,end));
            case 2
                fprintf('\nMean values of No. and ratio of sustained objects of all chunks:\t%d\t%0.2f%%\n\n', ...
                    ceil(mean(sustObjOfCurrChnkArr(1,:))),mean(sustObjOfCurrChnkArr(2,:)));
        end
    end

H.tElapsed = toc(SDCORstrt)-timWast; % Setting the elapsed time for SDCOR

guidata(hO,H);

end


%% Subfunctions goes Here!

function [bestSolParams,bestSolIdx,PSO_costArr] = PSO_DBSCAN(H)

global DBDS

lowBd = [0 2];
uppBd = [norm(max(DBDS)) H.dimCoef*H.p];

% % Setting the parameters bounds manually @-- debugging script --@
% lowBd = [1e2 1e2]; % @-- debugging script --@
% uppBd = [2e3 2e3]; % @-- debugging script --@

paramNo = numel(lowBd);
PSO_costArr = [];
particles = cell(0);
for c1 = 1:H.PSO_particleNo
    % Setting random values for epsilon and MinPts
    particles{1,c1} = unifrnd(lowBd,uppBd);
    % Calculating the cost of DBSCAN according to the parameters and
    % gaining the cluster indices
    [particles{2,c1},particles{3,c1}] = DBSCAN_Cost(particles{1,c1}(1),particles{1,c1}(2));
    % Setting the value for the velocity
    particles{4,c1} = unifrnd(-abs(uppBd-lowBd),abs(uppBd-lowBd));
end
% Setting the very first values for localBest
localBest = particles;

[minCost,minCostIdx] = min(cell2mat(particles(2,:)));
globalBest = particles(:,minCostIdx);

PSO_costArr = [PSO_costArr minCost];

H.PSO_finCond = 0;
plotOptional(H,{PSO_costArr},'PSOcost');

for c1 = 1:H.PSO_maxIter
    for c2 = 1:H.PSO_particleNo
        % Updating velocity
        particles{4,c2} = H.PSO_W*particles{4,c2}+ ...
            H.PSO_C1*rand(1,paramNo).*(localBest{1,c2}-particles{1,c2})+ ...
            H.PSO_C2*rand(1,paramNo).*(globalBest{1}-particles{1,c2});
        
        % Updating position
        particles{1,c2} = particles{1,c2}+particles{4,c2};
        
        % Checking boundaries
        uppBdCond = particles{1,c2}>uppBd;
        particles{1,c2}(uppBdCond) = uppBd(uppBdCond);
        lowBdCond = particles{1,c2}<lowBd;
        particles{1,c2}(lowBdCond) = lowBd(lowBdCond);
        
        % Updating the cost value
        [particles{2,c2},particles{3,c2}] = DBSCAN_Cost(particles{1,c2}(1),particles{1,c2}(2));
        
        % Updating the local and global solutions
        if particles{2,c2}<localBest{2,c2}
            localBest(:,c2) = particles(:,c2);
            if localBest{2,c2}<globalBest{2}
                globalBest = localBest(:,c2);
            end
        end
    end
    % Updating W
    H.PSO_W = (1-H.PSO_alpha)*H.PSO_W;
    
    PSO_costArr = [PSO_costArr globalBest{2}];
    
    plotOptional(H,{PSO_costArr},'PSOcost');
end

% Setting the best values for best found solution
globalBest{1}(2) = ceil(globalBest{1}(2));
bestSolParams = globalBest{1};
bestSolIdx = globalBest{3};

[costArrs] = costArrsCheck({PSO_costArr},1e100,10);
PSO_costArr = costArrs{1};

H.PSO_finCond = 1;
plotOptional(H,{PSO_costArr,bestSolParams},'PSOcost');

end

function [costArrs] = costArrsCheck(costArrs,maxVal,replcVal)

for c1 = 1:numel(costArrs)
    costArr = costArrs{c1};
    idxMaxVal = costArr==maxVal;
    if ~all(idxMaxVal)
        costArr(idxMaxVal) = max(costArr(~idxMaxVal))+replcVal;
    else
        costArr(idxMaxVal) = replcVal;
    end
    
    costArrs{c1} = costArr;
end

end

function [cost,idx] = DBSCAN_Cost(epsilon,MinPts)

[idx,~] = DBSCAN(epsilon,MinPts);

idxUnq = unique(idx(idx~=0));
K = numel(idxUnq);

if K>1
    %     cost = sum([DBindex(idx), CSindex(idx)]);
    cost = sum([DBindex(idx), BRindex(idx)]);
else
    cost = 1e100;
end

end

function [idx,isNoise] = DBSCAN(epsilon,MinPts)

global DBDS BLK_SZ_LIM

n = size(DBDS,1);

%------- error handing -------%
[~,sysView] = memory;
if n^2+BLK_SZ_LIM^2*8 > sysView.PhysicalMemory.Available
    % Creating the structure of the message box
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    msgCont = {'\fontsize{10}There is not {\bf{enough memory}} to run {\bf{DBSCAN}} on data available in RAM!, regarding the {\it{array-based implementation}} of MATLAB codes, and due to the {\it{incorrect input parameters}}, including random sampling rate {\bf{eta}}, or chunk size, or maybe DBSCAN parameters, {\bf{eps}} and {\bf{minPts}}.'; '\fontsize{10}Thus, please {\color{red}\bf{stop}} the execution and run it again with the true optimal parameters, or even just free up some memory space! ;-)'; '\newline{\color{red}\bf{Note:}} Consider the {\bf{BlckSzLim}} value in GUI as a parameter too!'};
    h = msgbox(msgCont,'Memory Error','error',CreateStruct);
    uiwait(h)
    keyboard % Stop the operation from here, because of lack of memory
end
%-----------------------------%

idx = sparse(n,1);
isNoise = sparse(n,1);
clustNo = 0;
neighbMat_Maker();
idxCore = find(sum(neighbMatX,2)>=MinPts);
corPntNo = numel(idxCore);
visited = sparse(n,1);

for d1 = 1:corPntNo
    if ~visited(idxCore(d1))
        clustNo = clustNo+1;
        clust = idxCore(d1);
        idxCorConsd = [];
        while true
            clustTemp = clust;
            clust = union(clust,find(any(neighbMatX(intersect(clust,setdiff(idxCore,idxCorConsd)),:),1)));
            if isempty(setdiff(clust,clustTemp))
                break;
            end
            idxCorConsd = intersect(idxCore,clustTemp);
        end
        idx(clust) = clustNo;
        visited(clust) = 1;
    end
end

isNoise(idx==0) = 1;
idx = full(idx);

% Nested function for computing the neighborhood matrix
    function neighbMat_Maker()
        
        neighbMatX = logical([]);
        
        cntLim = ceil(n/BLK_SZ_LIM);
        
        for c1 = 1:cntLim
            %     c1 % @-- debugging script --@
            if c1 ~= cntLim
                Y1 = DBDS((c1-1)*BLK_SZ_LIM+1:c1*BLK_SZ_LIM,:);
            else
                Y1 = DBDS((c1-1)*BLK_SZ_LIM+1:end,:);
            end
            neighbMatY = logical([]);
            
            for c2 = 1:cntLim
                %         c2 % @-- debugging script --@
                if c2 ~= cntLim
                    indVec = (c2-1)*BLK_SZ_LIM+1:c2*BLK_SZ_LIM;
                else
                    indVec = (c2-1)*BLK_SZ_LIM+1:n;
                end
                Y2 = DBDS(indVec,:);
                neighbMatY = [neighbMatY pdist2(Y1,Y2)<=epsilon];
            end
            
            neighbMatX = [neighbMatX; neighbMatY];
        end
        
    end

end

function [DBindex] = DBindex(idx)

idxUnq = unique(idx(idx~=0));
K = numel(idxUnq);

[cents,d2cMean] = centroidFind(idx);

centsDist = pdist2(cents,cents);

d2cMeanMat = repmat(d2cMean,1,K);

Ctemp = (d2cMeanMat+transpose(d2cMeanMat))./centsDist;
Ctemp(1:(K+1):K^2) = -inf;

DBindex = mean(max(Ctemp));

end

function [BRindex] = BRindex(idx)

[~,idxFreq,~] = idxFreqCalc(idx(idx~=0));

[~,d2cMean] = centroidFind(idx);

d2cMeanLog = log(d2cMean);
d2cMeanLog(isinf(d2cMeanLog)) = 0;
BRindex = sum(idxFreq.*d2cMeanLog);

end

function [cents,d2cMean] = centroidFind(idx)

global DBDS

dim = size(DBDS,2);
idxUnq = unique(idx(idx~=0));
K = numel(idxUnq);
cents = zeros(K,dim);
d2cMean = zeros(K,1);
for c1 = 1:K
    Y = DBDS(idx==idxUnq(c1),:);
    n1 = size(Y,1);
    cents(c1,:) = mean(Y);
    d2cMean(c1) = mean(sqrt(sum((Y-repmat(cents(c1,:),n1,1)).^2,2)));
end

end

function [idxUnq,idxFreq,K] = idxFreqCalc(idx)

if size(idx,1)==1
    idx = transpose(idx);
end
idxUnq = unique(idx);
K = length(idxUnq);
idxFreq = zeros(K,1);
for c1 = 1:K
    idxFreq(c1) = sum(idx==idxUnq(c1));
end

end

function [idxFin,singChck] = posSemiDefCheck(X,idx)

dim = size(X,2);
[idxUnq,idxFreq,K] = idxFreqCalc(idx(idx~=0));
idxFin = idx;

if K==0; singChck = 1; return; end % Check if no cluster is detected! (Regularly happens at 'Sampling' phase)

singChck = 0;
for c1 = 1:K
    clust = X(idx==idxUnq(c1),:);
    cCovEig = eig(cov(clust));
    
    cond1 = any(cCovEig<0);
    cond2 = ~isreal(cCovEig);
    cond3 = idxFreq(c1)<=dim;
    if cond1 || cond2 || cond3
        idxFin(idx==idxUnq(c1)) = 0;
        singChck = 1;
        
%         fprintf('Singularity happened!\n'); % @-- debugging script --@
        
    end
end

end

function clustInfMaker(H,data)

global CLUSTS

X = data{1};
idx = data{2};

currK = size(CLUSTS,2);
idxUnq = unique(idx(idx~=0));
K = numel(idxUnq);
for c1 = currK+1:currK+K
    Y = X(idx==idxUnq(c1-currK),:);
    n = size(Y,1);
    CLUSTS{1,c1} = mean(Y);
    CLUSTS{2,c1} = transpose(Y-repmat(CLUSTS{1,c1},n,1))*(Y-repmat(CLUSTS{1,c1},n,1));
    [coeff,~,latent,~,explained,~] = pca(Y);
    explCumSum = cumsum(explained); explCumSum(end) = 100;
    numComp = find(explCumSum>=H.PCvarRat*100,1);
    CLUSTS{3,c1} = coeff(:,1:numComp);
    CLUSTS{4,c1} = CLUSTS{1,c1}*CLUSTS{3,c1};
    CLUSTS{5,c1} = sqrt(transpose(latent(1:numComp)));
    CLUSTS{6,c1} = n;
    CLUSTS{7,c1} = numComp;
    CLUSTS{8,c1} = det((1/(n-1))*CLUSTS{2,c1});
end

end

function [clustBestArr] = clustAccMahal(H,X)

global CLUSTS CLUST_MODF_ARR CHNK_MEMB_COND

K = numel(CLUST_MODF_ARR);

mahalDistKarr = zeros(size(X,1),K);
for c1 = 1:K
    mahalDistKarr(:,c1) = sqrt(sum(((X*cell2mat(CLUSTS(3,CLUST_MODF_ARR(c1)))-cell2mat(CLUSTS(4,CLUST_MODF_ARR(c1))))./ ...
        cell2mat(CLUSTS(5,CLUST_MODF_ARR(c1)))).^2,2));
end

numCompArr = transpose(cell2mat(CLUSTS(7,CLUST_MODF_ARR)));

[mahalDists,clustBestArr] = min(mahalDistKarr,[],2);

clustBestArr(mahalDists>H.alphaMemb*sqrt(numCompArr(clustBestArr))) = 0;

if CHNK_MEMB_COND
    CLUST_MODF_ARR = unique(clustBestArr(clustBestArr~=0));
else
    clustBestArr(clustBestArr~=0) = CLUST_MODF_ARR(clustBestArr(clustBestArr~=0));
    CLUST_MODF_ARR = unique(clustBestArr(clustBestArr~=0));
end

end

function clustInfUpdate(X,clustBestArr,PCvarRat)

global CLUSTS

updClstK = unique(clustBestArr(clustBestArr~=0));
for c1 = 1:numel(updClstK)
    % Updating primary info. of temporary clusters
    Cmembs = X(clustBestArr==updClstK(c1),:);
    Csz = size(Cmembs,1);
    CLUSTS{2,updClstK(c1)} = CLUSTS{2,updClstK(c1)}+transpose(Cmembs-repmat(CLUSTS{1,updClstK(c1)},Csz,1))*...
        (Cmembs-repmat(CLUSTS{1,updClstK(c1)},Csz,1));
    CLUSTS{6,updClstK(c1)} = CLUSTS{6,updClstK(c1)}+Csz;
    
    % Updating secondary info. of temporary clusters
    [coeff,latent,explained] = pcacov((1/(CLUSTS{6,updClstK(c1)}-1))*CLUSTS{2,updClstK(c1)});
    explCumSum = cumsum(explained); explCumSum(end) = 100;
    CLUSTS{7,updClstK(c1)} = find(explCumSum>=PCvarRat*100,1);
    
    CLUSTS{3,updClstK(c1)} = coeff(:,1:CLUSTS{7,updClstK(c1)});
    CLUSTS{4,updClstK(c1)} = CLUSTS{1,updClstK(c1)}*CLUSTS{3,updClstK(c1)};
    CLUSTS{5,updClstK(c1)} = sqrt(transpose(latent(1:CLUSTS{7,updClstK(c1)})));
end

end

function [retainSet,retIdx] = retSetClustMembCheck(H,retainSet,retIdx)

global CLUST_MODF_ARR

if ~isempty(CLUST_MODF_ARR)
    while true
        [clustBestArrRet] = clustAccMahal(H,retainSet);
        
        clustInfUpdate(retainSet,clustBestArrRet,H.PCvarRat);
        
        idxRetIn = find(clustBestArrRet==0);
        
        retainSet = retainSet(idxRetIn,:);
        retIdx = retIdx(idxRetIn);
        
		if isempty(CLUST_MODF_ARR); break; end
    end
end

end

function [retainSet,retIdx] = retSetClustMaker(H,retainSet,retIdx)

global CLUSTS

retSz = size(retainSet,1);
if retSz~=0
    [idxRet,~] = DBSCAN(H.epsilonFin,H.MinPtsFin);
    [idxRetFin] = clustAccpDBSCAN();
    
    clustInfMaker(H,{retainSet,idxRetFin});
    
    retainSet = retainSet(idxRetFin==0,:);
    retIdx = retIdx(idxRetFin==0);
end

% Nested function for verifying DBSCAN output clusters
%     function [idxRetFin,unqFinLength] = clustAccpDBSCAN()
    function [idxRetFin] = clustAccpDBSCAN()
        
        [idxRet,~] = posSemiDefCheck(retainSet,idxRet);
        
        idxRetUnq = unique(idxRet(idxRet~=0));
        K = numel(idxRetUnq);
        
        idxRetFin = idxRet;
        for c1 = 1:K
            Z = retainSet(idxRet==idxRetUnq(c1),:);
            Zmean = mean(Z);
            Zdet = det(cov(Z));
            
            MDarr = zeros(1,H.origK);
            for c2 = 1:H.origK
                MDarr(c2) = sqrt(sum(((Zmean*cell2mat(CLUSTS(3,c2))-cell2mat(CLUSTS(4,c2)))./ ...
                    cell2mat(CLUSTS(5,c2))).^2));
            end
            [~,nrstInitClstIdx] = min(MDarr);
            
            
            if Zdet>H.sampDetArr(nrstInitClstIdx)
                
%                 fprintf('Delta Violation happened!\n'); % @-- debugging script --@
%                 beep % @-- debugging script --@
                
                [trueK,idxTrue] = findRetSetCorrK(Z,H.sampDetArr(nrstInitClstIdx),H.epsilonFin,H.MinPtsFin);
                
                if trueK~=1
                    idxRetFin(idxRetFin==idxRetUnq(c1)) = idxRetUnq(c1)+idxTrue/(trueK+1);
                    
%                     fprintf('Kmeans for irregular minicluster is utilized!\n'); % @-- debugging script --@
%                     beep % @-- debugging script --@
                    
                else
                    idxRetFin(idxRetFin==idxRetUnq(c1)) = 0;
                    
%                     fprintf('Kmeans for irregular minicluster failed!\n'); % @-- debugging script --@
%                     beep % @-- debugging script --@
                    
                end
            end
        end
        
    end

end

function [trueK,idxTrue] = findRetSetCorrK(X,deltaDet,epsilon,MinPts)

global DBDS

[n,p] = size(X);
idxArr = ones(n,1);
kMax = floor(n/(p+1));

trueK = [];
k = 2;
while true
    if k>kMax
        break;
    end
    
    idx = kmeans(X,k,'Replicates',5);
    idxArr = [idxArr idx];
    
    for c1 = 1:k
        PCDviol = det(cov(X(idx==c1,:)))>deltaDet;
        [~,singChck] = posSemiDefCheck(X(idx==c1,:),idx(idx==c1));
        
        DBDS = X(idx==c1,:);
        [idxKmns,~] = DBSCAN(epsilon,MinPts);
        DBSCANichr = numel(unique(idxKmns(idxKmns~=0)))>1;
        
        rejCond = PCDviol | singChck | DBSCANichr;
        
        if rejCond
            break;
        end
    end
    if ~rejCond
        trueK = k;
        break;
    else
        k = k+1;
        continue;
    end

end

if isempty(trueK)
    trueK = 1;
end
idxTrue = idxArr(:,trueK);

end

function [accResArr] = accReport(H,accResArr,retIdx)

retIdxFinal = sparse(H.n,1);
if ~isempty(retIdx)
    retIdxFinal(retIdx) = 1;
end

[~,~,~,AUC] = perfcurve(H.labFin,retIdxFinal,1);

TPno = sum(H.labFin(retIdx)==1);
FPno = sum(H.labFin(retIdx)==0);
FNno = H.OLno-TPno;
precision = TPno/(TPno+FPno);
recall = TPno/(TPno+FNno);
F1Meas = (2*precision*recall)/(precision+recall);

% fprintf('TPno = %d\tFPno = %d\tFNno = %d\n',TPno,FPno,FNno); % @-- debugging script --@

if isnan(precision); precision = 0; end
if isnan(recall); recall = 0; end
if isnan(F1Meas); F1Meas = 0; end

accResArr = [accResArr [AUC*100; precision*100; recall*100; F1Meas*100]];

plotOptional(H,{accResArr},'accPerChunk');

end

function [idxMeans,finalClusts,meansMeans,regenDS,idxRegenDS] = finalClustsMaker(H)

global CLUSTS

finalClusts = cell(0);

[idxMeans] = kmeans(H.means,H.origK,'Replicates',5);
[~,idxMeansFreq,~] = idxFreqCalc(idxMeans);

for c1 = 1:H.origK
    
    clustInf = CLUSTS([1 2 6],idxMeans==c1);
    
    if idxMeansFreq(c1)==1
        finalClusts{1,c1} = clustInf{1};
        n1 = clustInf{3};
        if n1<2
            unbiasCoeff = 1;
        else
            unbiasCoeff = 1/(n1-1);
        end
        finalClusts{2,c1} = unbiasCoeff*clustInf{2};
        sampNo = ceil(H.sampRate*n1);
        
        if sampNo<=H.p
            sampNo = H.p+1;
        end
        
		regenData = mvnrnd(clustInf{1},finalClusts{2,c1},sampNo);
		clustPrun();
        
        [coeff,latent,explained] = pcacov(finalClusts{2,c1});
        explCumSum = cumsum(explained); explCumSum(end) = 100;
        numComp = find(explCumSum>=H.PCvarRat*100,1);
        finalClusts{3,c1} = coeff(:,1:numComp);
        finalClusts{4,c1} = finalClusts{1,c1}*finalClusts{3,c1};
        finalClusts{5,c1} = sqrt(transpose(latent(1:numComp)));
        finalClusts{6,c1} = numComp;
        
        if H.dispOn
            finalClusts{7,c1} = [regenData c1*ones(size(regenData,1),1)];
        end
        
    else
        clustsSz = cell2mat(transpose(clustInf(3,:)));
        clustsMeans = cell2mat(transpose(clustInf(1,:)));
        finalClusts{1,c1} = sum(clustsSz.*clustsMeans)/sum(clustsSz);
        unbiasCoeffVec = 1./(clustsSz-1); unbiasCoeffVec(isinf(unbiasCoeffVec) | unbiasCoeffVec<0) = 1;
        sampNoVec = ceil(H.sampRate.*clustsSz); sampNoVec(sampNoVec<=H.p) = H.p+1;
        
        while true
            regenData = [];
            for c2 = 1:size(clustInf,2)
                clustCov = unbiasCoeffVec(c2)*clustInf{2,c2};
                regenData = [regenData; mvnrnd(clustInf{1,c2},clustCov,sampNoVec(c2))];
            end
            clustPrun();
            
            [~,singChck] = posSemiDefCheck(regenData,ones(size(regenData,1),1));
            if ~singChck
                break
            else
                fprintf('Singularity happened during regeneration process! The procedure will be conducted again\n');
                beep
            end
        end
        
        finalClusts{2,c1} = cov(regenData);
        
        [coeff,~,latent,~,explained,~] = pca(regenData);
        explCumSum = cumsum(explained); explCumSum(end) = 100;
        numComp = find(explCumSum>=H.PCvarRat*100,1);
        finalClusts{3,c1} = coeff(:,1:numComp);
        finalClusts{4,c1} = finalClusts{1,c1}*finalClusts{3,c1};
        finalClusts{5,c1} = sqrt(transpose(latent(1:numComp)));
        finalClusts{6,c1} = numComp;
        
        if H.dispOn
            finalClusts{7,c1} = [regenData c1*ones(size(regenData,1),1)];
        end
        
    end
    
end

if H.dispOn
    meansMeans = cell2mat(transpose(finalClusts(1,:)));
    regenDS = cell2mat(transpose(finalClusts(7,:)));
    idxRegenDS = regenDS(:,end);
    regenDS = regenDS(:,1:end-1);
else
    meansMeans = [];
    regenDS = [];
    idxRegenDS = [];
end

finalClusts = finalClusts(1:6,:);


%% Nested functions here!

% Nested function for pruning the regenerated data
    function clustPrun()
        mahalDist = sqrt(mahal(regenData,regenData));
        regenData = regenData(mahalDist<=H.betaPrun*sqrt(H.p),:);
    end

end

function [mahalScores,idxFin,finalAUC] = OLscoreMaker(H)

K = size(H.finalClusts,2);
mahalScorKarr = zeros(H.n,K);
maxIter = ceil(H.n/H.chunkSz);

for c1 = 1:maxIter
    indStart = (c1-1)*H.chunkSz+1;
    if c1~=maxIter
        indEnd = c1*H.chunkSz;
    else
        indEnd = H.n;
    end
    
    for c2 = 1:K
        mahalScorKarr(indStart:indEnd,c2) = sqrt(sum(((H.DS(indStart:indEnd,:)*cell2mat(H.finalClusts(3,c2))-cell2mat(H.finalClusts(4,c2)))./ ...
            cell2mat(H.finalClusts(5,c2))).^2,2));
    end
end

[mahalScores,idxFin] = min(mahalScorKarr,[],2);

[~,~,~,finalAUC] = perfcurve(H.labFin,mahalScores,1);

set(H.finalAUCbyScores_statText,'String',num2str(finalAUC,'%0.3f'));

end

