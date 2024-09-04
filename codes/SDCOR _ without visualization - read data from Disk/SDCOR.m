
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ce.aut.ac.ir/~sann_cv/
% June 2020

function [] = SDCOR(hObject,handles)

    SDCORstrt = tic; % Execution time for SDCOR
    
    global DBDS CLUSTS CLUST_MODF_ARR CHNK_MEMB_COND
    PSO_time = 0;
    
    %------- displaying progress -------%
    handles.progLevl_statText.FontSize = 13;
    handles.progLevl_statText.String = 'Sampling phase started!'; pause(.001);
    %-----------------------------------%
    
    % tSt = tic; % Setting start time for sampling phase @-- debugging script --@
    rndSmp_Maker();
    while true
        sampPrmChosMthd();
        [handles.idxSamp,singChck] = posDefCheck(handles.sampData,handles.idxSamp);
        
        if ~singChck
            break;
        else
            fprintf('Singularity happened during "Sampling" phase! The procedure requires to be conducted again.\n'); beep
            keyboard
        end
    end
    % tElp0 = toc(tSt); % Setting elapsed time for sampling phase @-- debugging script --@
    % fprintf('\n\nRandSamp:\t%0.2fe-4 sec\n\n',tElp0*1e4); % @-- debugging script --@
    
    % Initializing the temporary clustering model
    CLUSTS = cell(0);
    clustInfMaker(handles,{handles.sampData,handles.idxSamp,handles.idxSamp});
    handles.origK = size(CLUSTS,2);
    handles.origK_val_statText.String = num2str(handles.origK);
    
    % Building the initial determinant array
    handles.sampDetArr = zeros(1,handles.origK);
    for d1 = 1:handles.origK
        handles.sampDetArr(d1) = det((1/(CLUSTS{6,d1}-1))*CLUSTS{2,d1});
    end
    handles.sampDetArr(handles.sampDetArr<=0) = 1; % Error handling
    
    handles = rmfield(handles,{'sampData'}); % Clearing RAM from unnecessary data
    
    %------- displaying progress -------%
    handles.progLevl_statText.String = 'Sampling phase finished!'; pause(.5);
    %-----------------------------------%
    
    % % Making the array of the number and ratio of sustained objects in RAM belonging to each chunk @-- debugging script --@
    % sustObjOfCurrChnkArr = []; % @-- debugging script --@
    
    retainSet = [];
    retIdx = [];
    maxIter = ceil(handles.n/handles.chunkSz);
    
    %------- displaying progress -------%
    handles.progLevl_statText.String = '';
    handles.progLevl_statText.FontSize = 15;
    handles.progLevl_statText.String = [num2str(0) '/' num2str(maxIter)]; pause(.001);
    %-----------------------------------%
    
    for c1 = 1:maxIter
        indStart = (c1-1)*handles.chunkSz+1;
        if c1~=maxIter
            indEnd = c1*handles.chunkSz;
        else
            indEnd = handles.n;
        end
        chunk = handles.labDS.X(indStart:indEnd,:);
        
    %     tSt = tic; % Setting start time for chunk memb. @-- debugging script --@
        CHNK_MEMB_COND = 1;
        CLUST_MODF_ARR = transpose(1:size(CLUSTS,2));
        [clustBestArr] = clustAccMahal(handles,chunk);
        CHNK_MEMB_COND = 0;
        
        clustInfUpdate(chunk,clustBestArr,handles.PCvarRat);
    %     tElp1 = toc(tSt); % Setting elapsed time for chunk memb. @-- debugging script --@
        
        retainSet = [retainSet; chunk(clustBestArr==0,:)];
        retIdx = [retIdx; (c1-1)*handles.chunkSz+find(clustBestArr==0)];
        
    %     tElp2 = 0; tElp3 = 0; tElp4 = 0; % @-- debugging script --@
        if ~isempty(retIdx)
    %         tSt = tic; % Setting start time for retain set memb. 1 @-- debugging script --@
            [retainSet,retIdx] = retSetClustMembCheck(handles,retainSet,retIdx);
    %         tElp2 = toc(tSt); % Setting elalsed time for retain set memb. 1 @-- debugging script --@
            
		    if ~isempty(retIdx)
    % 	        tSt = tic; % Setting start time for retain set clust. @-- debugging script --@
			    befrLastCretClstNoK = size(CLUSTS,2);
			    DBDS = retainSet;
			    [retainSet,retIdx] = retSetClustMaker(handles,retainSet,retIdx);
                CLUST_MODF_ARR = transpose(befrLastCretClstNoK+1:size(CLUSTS,2));
    % 	        tElp3 = toc(tSt); % Setting elapsed time for retain set clust. @-- debugging script --@
			    
                if ~isempty(retIdx)
    %                 tSt = tic; % Setting start time for retain set memb. 2 @-- debugging script --@
                    [retainSet,retIdx] = retSetClustMembCheck(handles,retainSet,retIdx);
    %                 tElp4 = toc(tSt); % Setting elalsed time for retain set memb. 2 @-- debugging script --@
                end
		    end
            
        end
	    
    %     fprintf('ChunkMemb:\t%0.2fe-4 sec;\tretSetMemb.1:\t%0.2fe-4 sec;\tretSetClust.:\t%0.2fe-4 sec;\tretSetMemb.2:\t%0.2fe-4 sec\n',...
    %         tElp1*1e4,tElp2*1e4,tElp3*1e4,tElp4*1e4); % Printing elapsed times for various stages @-- debugging script --@
        
    %     % Demonstrating the number and ratio of sustained objects in RAM belonging to current chunk @-- debugging script --@
    %     sustObjNoRAM(1); % @-- debugging script --@
        
        %------- displaying progress -------%
        handles.progLevl_statText.String = [num2str(c1) '/' num2str(maxIter)]; pause(.001);
        %-----------------------------------%
    end
    
    clearvars chunk DBDS retainSet retIdx clustBestArr % Clearing RAM from unnecessary data
    
    % % Demonstrating the mean value of numbers and ratios of sustained objects in RAM belonging to all processed chunks @-- debugging script --@
    % sustObjNoRAM(2); % @-- debugging script --@
    
    %------- displaying progress -------%
    handles.progLevl_statText.String = '';
    handles.progLevl_statText.FontSize = 10;
    handles.progLevl_statText.String = 'Building the final clustering model!'; pause(.5);
    %-----------------------------------%
    
    % tSt = tic; % Setting start time for building the final clustering model @-- debugging script --@
    [handles.finalClusts] = finalClustsMaker(handles);
    % tElp5 = toc(tSt); % Setting elapsed time for building the final clustering model @-- debugging script --@
    % fprintf('FinClstMakr:\t%0.2fe-4 sec\n',tElp5*1e4); % @-- debugging script --@
    
    clearvars CLUSTS % Clearing RAM from unnecessary data
    
    %------- displaying progress -------%
    handles.progLevl_statText.String = 'Calculating the outlierness scores!'; pause(.001);
    %-----------------------------------%
    
    % tSt = tic; % Setting start time for the scoring phase @-- debugging script --@
    [handles.mahalScores,handles.idxFin,handles.ROC,handles.PR] = OLscoreMaker(handles);
    % tElp6 = toc(tSt); % Setting elapsed time for the scoring phase @-- debugging script --@
    % fprintf('OLscrMakr:\t%0.2fe-4 sec\n',tElp6*1e4); % @-- debugging script --@
    
    %------- displaying progress -------%
    handles.progLevl_statText.String = '';
    handles.progLevl_statText.FontSize = 15;
    handles.progLevl_statText.String = 'D O N E !'; pause(.001);
    %-----------------------------------%
    
    handles = rmfield(handles,{'finalClusts'}); % Clearing RAM from unnecessary data
    
    %% Nested functions here!
    
    % Nested function for random sampling
        function rndSmp_Maker()
            if ~handles.nonUnifSamp
                % fix-step sampling with variation on the starting point. matfile()
                % function does not support randperm() indexing
                stepLeng = floor(1/handles.sampRate);
                handles.sampInd = transpose(randi(stepLeng):stepLeng:handles.n);
                handles.sampData = handles.labDS.X(handles.sampInd,:);
            
            else
                MatFileNonUnifSamp(hObject, handles, 'X');
                handles = guidata(hObject); % Retrieve the updated handles structure
    
            end
            
        end
    
    % Nested function for selecting the type of DBSCAN parameter choosing
        function sampPrmChosMthd()
            switch handles.PCM
                case 'Kgraph_pcm_radioBtn'
                    [~,KdistGrph] = knnsearch(handles.sampData,handles.sampData,'K',handles.manuMnPt);
                    KdistGrph = sort(KdistGrph(:,end),'descend');
                    
                    % plotting sorted k-dist graph
                    axes(handles.axes1); stairs(KdistGrph); legend('{\it{k}}-dist graph');
                    
                    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
                    msgCont = ['\fontsize{10} After closing this dialog box, the program will be suspended by the "keyboard" command. '...
                        'Please cease the operation by pressing "Shift+F5" and then return to the program window. Thereafter, first choose the '...
                        '"Data Cursor" tool from the "Axes Toolbar" above the window, and then in the "DBSCAN Param Choosing" axes, '...
                        'select a position which is indicating the first point in the first {\bf{"valley"}} of the sorted '...
                        '{\it{k}}-dist graph. Remember the Y value of this bar as the DBSCAN Epsilon parameter for clustering the sampled data. '...
                        '\newline\newline Finally, reset the whole disabled environment by pressing the "RESET BTNS" button; after that, '...
                        'set the "DBSCAN parameter choosing" approach in the program window as "Manual" and use the obtained quantity '...
                        'out of the {\it{k}}-dist graph for the sampling Epsilon parameter.\newline\newline {\color{red}\bf{Note:}} '...
                        'For the plots related to large datasets, you might need to zoom in to better find the first valley. Moreover, '...
                        'to evade the singularity occurrence, you will need to select a higher threshold through the same process '...
                        'already mentioned; even it might be more reasonable to start from lower values until you reach the minimum '...
                        'acceptable Epsilon (the best case).'];
                    h = msgbox(msgCont,'Sampling Epsilon selection','help',CreateStruct); uiwait(h)
                    keyboard % Stop the operation from here for visually selecting the best sampling Epsilon value
                    
                case 'PSO_pcm_radioBtn'
                    DBDS = handles.sampData;
                    [handles.paramSampDS,handles.idxSamp,PSO_costArr,PSO_time] = PSO_DBSCAN(handles);
                    
                    handles.paramCostArrSamp = {PSO_costArr};
                    handles.origEps = handles.epsCoeff*handles.paramSampDS(1);
                    handles.origMnPt = handles.paramSampDS(2);
                    
                case 'manu_pcm_radioBtn'
                    DBDS = handles.sampData;
                    [handles.idxSamp,~] = DBSCAN(handles.manuEps,handles.manuMnPt);
                    
                    handles.paramSampDS = [handles.manuEps handles.manuMnPt];
                    handles.paramCostArrSamp = {}; PSO_time = 0;
                    
                    handles.origEps = handles.epsCoeff*handles.paramSampDS(1);
                    handles.origMnPt = handles.paramSampDS(2);
            end
        end
    
    % Nested function for displaying the number of sustained objects of current chunk
        function sustObjNoRAM(iterType)
            switch iterType
                case 1
                    sustObjOfCurrChnk = sum(clustBestArr==0);
                    sustObjOfCurrChnkArr = [sustObjOfCurrChnkArr [sustObjOfCurrChnk; sustObjOfCurrChnk/handles.chunkSz*100]];
                    fprintf('No. and ratio of sustained objects of chunk %d is:\t%d\t%0.2f%%\n',c1,sustObjOfCurrChnkArr(1,end), ...
                        sustObjOfCurrChnkArr(2,end));
                case 2
                    fprintf('\nMean values of No. and ratio of sustained objects of all chunks:\t%d\t%0.2f%%\n\n', ...
                        ceil(mean(sustObjOfCurrChnkArr(1,:))),mean(sustObjOfCurrChnkArr(2,:)));
            end
        end
    
    handles.tElapsed = toc(SDCORstrt)-PSO_time-1-(maxIter+4)*0.001; % Setting the elapsed time for SDCOR
    
    guidata(hObject,handles);

end


%% Subfunctions goes Here!

function MatFileNonUnifSamp(hObject, handles, varName)
    % Get the variable size out of the MAT-file object
    matObj = handles.labDS;
    varSize = size(matObj, varName);

    % Calculate the available memory in MATLAB (in bytes)
    userview = memory; 
    maxRAM = userview.MemAvailableAllArrays;
    divisionsNo = 2;
    maxRAM = maxRAM / divisionsNo; % Consuming a portion of the available memory to avoid probable further conflicts

    % Calculate the maximum number of rows that fit in RAM
    bytesPerUnit = 8; % Size of the numeric data type (class) in bytes. Here we use double which uses 8 bytes.
    rowSize = varSize(2) * bytesPerUnit;  % Each double is 8 bytes
    maxRowsThatFitInRAM = floor(maxRAM / rowSize);

    % Initialize an empty array to collect sampled data and global indices
    handles.sampData = [];
    handles.sampInd = [];

    % Loop through the data in chunks
    for startRow = 1:maxRowsThatFitInRAM:varSize(1)
        endRow = min(startRow + maxRowsThatFitInRAM - 1, varSize(1));

        % Load the chunk of data into memory
        partialData = matObj.(varName)(startRow:endRow, :);

        % Perform non-uniform sampling based on the number of rows in the chunk
        [sampledInd, sampIndGlobal] = nonUniformSampling(size(partialData, 1), handles.sampRate, startRow);

        % Collect the sampled data from this chunk
        sampledChunk = partialData(sampledInd, :);
        handles.sampData = [handles.sampData; sampledChunk];  % Append sampled data from this chunk

        % Collect global sampling indices
        handles.sampInd = [handles.sampInd; sampIndGlobal];

        % Clear partialData from memory
        clear partialData;
    end

    guidata(hObject,handles);

end

function [sampledInd, sampIndGlobal] = nonUniformSampling(numRows, sampRate, startRow)
    % Non-uniform sampling across multiple sections
    sectNo = 5;  % Number of sections to divide the dataset for non-uniform sampling
    randSect = rand(1, sectNo);  % Generate random values for each section
    randSect = randSect / sum(randSect);  % Normalize the random values so their sum equals 1
    randSect = randSect * 0.75 + 0.25 / sectNo;  % Ensure a minimum sampling size by adjusting the section sizes

    % Calculate the starting indices for each section based on cumulative sum
    sampStartInd = floor(cumsum(randSect) * numRows);
    sampStartInd = [0, sampStartInd(1:end-1) - 1];  % Adjust starting indices for indexing

    % Sort sections by size and flip for sampling share allocation
    [randSectSort, randSectSortIdx] = sort(randSect);  % Sort sections by size
    randSectShare = flip(randSectSort);  % Flip the sorted sections to give the largest section the smallest share
    randSectShare(randSectSortIdx) = randSectShare;  % Reassign the shares to their original sections

    % Determine the number of samples for each section
    randSectSize = floor(randSect * numRows);  % Calculate the number of samples in each section
    randSectShareSize = floor(randSectShare * sampRate * numRows);  % Calculate the sampling share size for each section

    sampledInd = [];  % Initialize an empty array to collect sampling indices from all sections

    % Perform random sampling within each section
    for sec = 1:sectNo
        sampledInd = [sampledInd, randperm(randSectSize(sec), randSectShareSize(sec)) + sampStartInd(sec)];
    end

    % Convert the sampling indices to a column vector
    sampledInd = sampledInd';
    
    % Calculate global indices corresponding to the original dataset
    sampIndGlobal = sampledInd + (startRow - 1);

end

function [bestSolParams,bestSolIdx,PSO_costArr,PSO_time] = PSO_DBSCAN(handles)

    tic
    global DBDS
    
    minMnPt = floor(log(handles.n));
    maxMnPt = handles.manuMnPt; if maxMnPt<minMnPt; maxMnPt = minMnPt; end
    [~,Kdist] = knnsearch(DBDS,DBDS,'K',minMnPt); Kdist = Kdist(:,end);
    lowBnd = [min(Kdist) minMnPt];
    uppBnd = [max(Kdist) maxMnPt];
    
    paramNo = numel(lowBnd);
    PSO_costArr = [];
    particles = cell(0);
    for c1 = 1:handles.PSO_particleNo
        % Setting random values for epsilon and MinPts
        particles{1,c1} = unifrnd(lowBnd,uppBnd);
        % Calculating the cost of DBSCAN according to the parameters and
        % gaining the cluster indices
        [particles{2,c1},particles{3,c1}] = DBSCAN_Cost(particles{1,c1}(1),ceil(particles{1,c1}(2)));
        % Setting a random value for the velocity
        particles{4,c1} = unifrnd(-abs(uppBnd-lowBnd),abs(uppBnd-lowBnd));
    end
    % Setting the very first values for localBest
    localBest = particles;
    
    [minCost,minCostIdx] = min(cell2mat(particles(2,:)));
    globalBest = particles(:,minCostIdx);
    
    PSO_costArr = [PSO_costArr minCost];
    
    handles.PSO_finCond = 0;
    PSOplot({PSO_costArr});
    
    for c1 = 1:handles.PSO_maxIter
        for c2 = 1:handles.PSO_particleNo
            % Updating velocity
            particles{4,c2} = handles.PSO_W*particles{4,c2}+ ...
                handles.PSO_C1*rand(1,paramNo).*(localBest{1,c2}-particles{1,c2})+ ...
                handles.PSO_C2*rand(1,paramNo).*(globalBest{1}-particles{1,c2});
            
            % Updating position
            particles{1,c2} = particles{1,c2}+particles{4,c2};
            
            % Checking boundaries
            uppBndCond = particles{1,c2}>uppBnd; particles{1,c2}(uppBndCond) = uppBnd(uppBndCond);
            lowBndCond = particles{1,c2}<lowBnd; particles{1,c2}(lowBndCond) = lowBnd(lowBndCond);
            
            
            % Updating the cost value
            [particles{2,c2},particles{3,c2}] = DBSCAN_Cost(particles{1,c2}(1),ceil(particles{1,c2}(2)));
            
            % Updating the local and global solutions
            if particles{2,c2}<localBest{2,c2}
                localBest(:,c2) = particles(:,c2);
                if localBest{2,c2}<globalBest{2}
                    globalBest = localBest(:,c2);
                end
            end
        end
        % Updating W
        handles.PSO_W = (1-handles.PSO_alpha)*handles.PSO_W;
        
        % Updating and plotting the PSO cost array
        PSO_costArr = [PSO_costArr globalBest{2}];
        PSOplot({PSO_costArr});
    end
    
    % Setting the best values for the best found solution
    globalBest{1}(2) = ceil(globalBest{1}(2));
    bestSolParams = globalBest{1};
    bestSolIdx = globalBest{3};
    
    % Setting inf values as a higher value than the maximum for better plotting
    PSO_costArr(isinf(PSO_costArr)) = 10*max(PSO_costArr);
    
    handles.PSO_finCond = 1;
    PSOplot({PSO_costArr,bestSolParams});
    
    PSO_time = toc; % computing the wasted time for PSO
    
    % Nested function for plotting PSO cost array
        function PSOplot(data)
            axes(handles.axes1);
            plot(data{1},'-r');
            pause(.001);
            if handles.PSO_finCond
                legend(sprintf('PSO costArr for SampDS\nEps=%0.3f, MinPts=%d',data{2}(1),data{2}(2)),'location','best'); grid on;
            end
        end

end

function [cost,idx] = DBSCAN_Cost(epsilon,MinPts)

    global DBDS
    
    [idx,~] = DBSCAN(epsilon,MinPts);
    
    [idx,singChck] = posDefCheck(DBDS,idx);
    if ~singChck && any(idx==0)
        cost = sum([DBindex(idx), CSindex(idx), sum(idx==0)/size(DBDS,1)]);
    else
        cost = inf;
    end

end

function [idx,corepts] = DBSCAN(epsilon,MinPts)

    global DBDS
    [idx,corepts] = dbscan(DBDS,epsilon,MinPts);
    idx(idx==-1) = 0;  % We handle noise here with index 0 in lieu of typical -1

end

function [DBidx] = DBindex(idx)

    idxUnq = unique(idx);
    K = numel(idxUnq);
    
    [cents,d2cMean] = centroidFind(idx);
    centsDist = pdist2(cents,cents);
    d2cMeanMat = repmat(d2cMean,1,K);
    
    Ctemp = (d2cMeanMat+transpose(d2cMeanMat))./centsDist;
    Ctemp(1:(K+1):K^2) = -inf;
    DBidx = mean(max(Ctemp));

end

function [CSidx] = CSindex(idx)

    global DBDS
    
    idxUnq = unique(idx);
    K = numel(idxUnq);
    
    distMaxMean = zeros(K,1);
    for c1 = 1:K
        idxKgi = idx==idxUnq(c1);
        distMaxK = max(pdist2(DBDS(idxKgi,:),DBDS(idxKgi,:)),[],2);
        distMaxMean(c1) = mean(distMaxK);
    end
    
    [cents,~] = centroidFind(idx);
    centsDist = pdist2(cents,cents);
    centsDist(1:(K+1):K^2) = inf;
    centsDistMin = min(centsDist,[],2);
    
    CSidx = sum(distMaxMean)/sum(centsDistMin);

end

function [cents,d2cMean] = centroidFind(idx)

    global DBDS
    
    dim = size(DBDS,2);
    idxUnq = unique(idx);
    K = numel(idxUnq);
    cents = zeros(K,dim);
    d2cMean = zeros(K,1);
    for c1 = 1:K
        Y = DBDS(idx==idxUnq(c1),:);
        n = size(Y,1);
        cents(c1,:) = mean(Y,1);
        d2cMean(c1) = mean(sqrt(sum((Y-repmat(cents(c1,:),n,1)).^2,2)));
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

function [idxFin,singChck] = posDefCheck(X,idx)

    dim = size(X,2);
    [idxUnq,idxFreq,K] = idxFreqCalc(idx(idx~=0));
    idxFin = idx;
    
    if K==0; singChck = 1; return; end % Check if no cluster is detected! (Regularly happens at 'Sampling' phase)
    
    singChck = 0;
    for c1 = 1:K
        clust = X(idx==idxUnq(c1),:);
        cCovEig = eig(cov(clust));
        
        cond1 = any(cCovEig<=0);
        cond2 = ~isreal(cCovEig);
        cond3 = idxFreq(c1)<=dim;
        if cond1 || cond2 || cond3
            idxFin(idx==idxUnq(c1)) = 0;
            singChck = 1;
            
    %         fprintf('Singularity happened!\n'); % @-- debugging script --@
            
        end
    end

end

function clustInfMaker(handles,data)

    global CLUSTS
    
    X = data{1};
    idx = data{2};
    Kind = data{3};
    
    currK = size(CLUSTS,2);
    idxUnq = unique(idx(idx~=0));
    K = numel(idxUnq);
    for c1 = currK+1:currK+K
        idxKgi = idx==idxUnq(c1-currK);
        Y = X(idxKgi,:);
        n = size(Y,1);
        CLUSTS{1,c1} = mean(Y,1);
        CLUSTS{2,c1} = transpose(Y-repmat(CLUSTS{1,c1},n,1))*(Y-repmat(CLUSTS{1,c1},n,1));
        [coeff,~,latent,~,explained,~] = pca(Y);
        explCumSum = cumsum(explained); explCumSum(end) = 100;
        numComp = find(explCumSum>=handles.PCvarRat*100,1);
        CLUSTS{3,c1} = coeff(:,1:numComp);
        CLUSTS{4,c1} = CLUSTS{1,c1}*CLUSTS{3,c1};
        CLUSTS{5,c1} = sqrt(transpose(latent(1:numComp)));
        CLUSTS{6,c1} = n;
        CLUSTS{7,c1} = numComp;
        CLUSTS{8,c1} = unique(Kind(idxKgi));
    end

end

function [clustBestArr] = clustAccMahal(handles,X)

    global CLUSTS CLUST_MODF_ARR CHNK_MEMB_COND
    
    K = numel(CLUST_MODF_ARR);
    
    mahalDistKarr = zeros(size(X,1),K);
    for c1 = 1:K
        mahalDistKarr(:,c1) = sqrt(sum(((X*CLUSTS{3,CLUST_MODF_ARR(c1)}-CLUSTS{4,CLUST_MODF_ARR(c1)})./ ...
            CLUSTS{5,CLUST_MODF_ARR(c1)}).^2,2));
    end
    
    numCompArr = transpose(cell2mat(CLUSTS(7,CLUST_MODF_ARR)));
    
    [mahalDists,clustBestArr] = min(mahalDistKarr,[],2);
    
    clustBestArr(mahalDists>handles.alphaMemb*sqrt(numCompArr(clustBestArr))) = 0;
    
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

function [retainSet,retIdx] = retSetClustMembCheck(handles,retainSet,retIdx)

    global CLUST_MODF_ARR
    
    if ~isempty(CLUST_MODF_ARR)
        while true
            [clustBestArrRet] = clustAccMahal(handles,retainSet);
            
            clustInfUpdate(retainSet,clustBestArrRet,handles.PCvarRat);
            
            idxRetIn = find(clustBestArrRet==0);
            
            retainSet = retainSet(idxRetIn,:);
            retIdx = retIdx(idxRetIn);
            
		    if isempty(CLUST_MODF_ARR); break; end
        end
    end

end

function [retainSet,retIdx] = retSetClustMaker(handles,retainSet,retIdx)

    global CLUSTS
    
    retSz = size(retainSet,1);
    if retSz~=0
        [idxRet,~] = DBSCAN(handles.origEps,handles.origMnPt);
        Kind = zeros(numel(idxRet),1);
        clustAccpDBSCAN();
        
        clustInfMaker(handles,{retainSet,idxRet,Kind});
        
        retainSet = retainSet(idxRet==0,:);
        retIdx = retIdx(idxRet==0);
    end
    
    % Nested function for verifying DBSCAN output clusters
        function clustAccpDBSCAN()
            
            [idxRet,~] = posDefCheck(retainSet,idxRet);
            
            idxRetUnq = unique(idxRet(idxRet~=0));
            K = numel(idxRetUnq);
            for c1 = 1:K
                idxRetKgi = idxRet==idxRetUnq(c1);
                Z = retainSet(idxRetKgi,:);
                Zmean = mean(Z);
                Zdet = det(cov(Z));
                
                MDarr = zeros(1,handles.origK);
                for c2 = 1:handles.origK
                    MDarr(c2) = sqrt(sum(((Zmean*CLUSTS{3,c2}-CLUSTS{4,c2})./ ...
                        CLUSTS{5,c2}).^2));
                end
                [~,nrstInitClstIdx] = min(MDarr);
                Kind(idxRetKgi) = nrstInitClstIdx;
                
                if Zdet>handles.sampDetArr(nrstInitClstIdx)
                    
    %                 fprintf('Delta Violation happened!\n'); % @-- debugging script --@
    %                 beep % @-- debugging script --@
                    
                    [trueK,idxTrue] = findRetSetCorrK(Z,handles.sampDetArr(nrstInitClstIdx),handles.origEps,handles.origMnPt);
                    
                    if trueK~=1
                        idxRet(idxRetKgi) = idxRetUnq(c1)+idxTrue/(trueK+1);
                        
    %                     fprintf('Kmeans for irregular minicluster is utilized!\n'); % @-- debugging script --@
    %                     beep % @-- debugging script --@
                        
                    else
                        idxRet(idxRetKgi) = 0;
                        Kind(idxRetKgi) = 0;
                        
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
            DBDS = X(idx==c1,:);
            [idxKmns,~] = DBSCAN(epsilon,MinPts);
            DBSCANichr = numel(unique(idxKmns(idxKmns~=0)))>1;
            
		    [~,singChck] = posDefCheck(X(idx==c1,:),idx(idx==c1));
            
		    PCDviol = det(cov(X(idx==c1,:)))>deltaDet;
            
            rejCond = DBSCANichr | singChck | PCDviol;
            
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

function [finalClusts] = finalClustsMaker(handles)

    global CLUSTS
    
    finalClusts = cell(0);
    idxMeans = cell2mat(CLUSTS(8,:));
    [~,idxMeansFreq,~] = idxFreqCalc(idxMeans);
    
    for c1 = 1:handles.origK
        
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
            
            [coeff,latent,explained] = pcacov(finalClusts{2,c1});
            explCumSum = cumsum(explained); explCumSum(end) = 100;
            numComp = find(explCumSum>=handles.PCvarRat*100,1);
            finalClusts{3,c1} = coeff(:,1:numComp);
            finalClusts{4,c1} = finalClusts{1,c1}*finalClusts{3,c1};
            finalClusts{5,c1} = sqrt(transpose(latent(1:numComp)));
            finalClusts{6,c1} = numComp;
            
        else
            clustsSz = cell2mat(transpose(clustInf(3,:)));
            clustsMeans = cell2mat(transpose(clustInf(1,:)));
            finalClusts{1,c1} = sum(clustsSz.*clustsMeans)/sum(clustsSz);
            unbiasCoeffVec = 1./(clustsSz-1); unbiasCoeffVec(isinf(unbiasCoeffVec) | unbiasCoeffVec<0) = 1;
            sampNoVec = ceil(handles.sampRate.*clustsSz); sampNoVec(sampNoVec<=handles.p) = handles.p+1;
            
            while true
                regenData = [];
                for c2 = 1:size(clustInf,2)
                    clustCov = unbiasCoeffVec(c2)*clustInf{2,c2};                
                    regenData = [regenData; mvnrnd(clustInf{1,c2},clustCov,sampNoVec(c2))];
                end
                
                % pruning the regenerated data
                mahalDist = sqrt(mahal(regenData,regenData));
                regenData = regenData(mahalDist<=handles.betaPrun*sqrt(handles.p),:);
                
                [~,singChck] = posDefCheck(regenData,ones(size(regenData,1),1));
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
            numComp = find(explCumSum>=handles.PCvarRat*100,1);
            finalClusts{3,c1} = coeff(:,1:numComp);
            finalClusts{4,c1} = finalClusts{1,c1}*finalClusts{3,c1};
            finalClusts{5,c1} = sqrt(transpose(latent(1:numComp)));
            finalClusts{6,c1} = numComp;
            
        end
        
    end
    
    finalClusts = finalClusts(1:6,:);

end

function [mahalScores,idxFin,ROC,PR] = OLscoreMaker(handles)

    K = size(handles.finalClusts,2);
    mahalScorKarr = zeros(handles.n,K);
    maxIter = ceil(handles.n/handles.chunkSz);
    
    for c1 = 1:maxIter
        indStart = (c1-1)*handles.chunkSz+1;
        if c1~=maxIter
            indEnd = c1*handles.chunkSz;
        else
            indEnd = handles.n;
        end
        
        for c2 = 1:K
            mahalScorKarr(indStart:indEnd,c2) = sqrt(sum(((handles.labDS.X(indStart:indEnd,:)*handles.finalClusts{3,c2}-handles.finalClusts{4,c2})./...
                handles.finalClusts{5,c2}).^2,2));
        end
    end
    
    [mahalScores,idxFin] = min(mahalScorKarr,[],2);
    
    [~,~,~,ROC] = perfcurve(handles.labDS.y,mahalScores,1);
    [~,~,~,PR] = perfcurve(handles.labDS.y,mahalScores,1,'XCrit','reca','YCrit','prec');

end

