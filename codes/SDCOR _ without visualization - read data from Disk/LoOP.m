
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ceit.aut.ac.ir/~sann_cv/
% June 2020

function [LoOPvals,AUC_LoOP,tElapsed] = LoOP(H)

global BLK_SZ_LIM
kVal = H.minPtsIntv(1);
lambda = H.LoOPlambda;

tStart = tic; % Set start time

matlMat_Maker();

kDistX = matlMatX(:,kVal);
kDistNeighbLog = matlMatX <= kDistX;
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
LoOPvals(LoOPvals < 0) = 0;

tElapsed = toc(tStart); % Set end time

% Calculating AUC
[~,~,~,AUC_LoOP] = perfcurve(H.labDS.y,LoOPvals,1);


%% Nested Functions Here!

% Nested function for calculating the neighborhood graph block-by-block
    function matlMat_Maker()
        
        [n,~] = size(H.labDS,'X');
        
        %------- Error handing -------%
        [~,sysView] = memory;
        if n*kVal*8+BLK_SZ_LIM^2*8*6+BLK_SZ_LIM*kVal*8*3 > sysView.PhysicalMemory.Available
            % Creating the structure of the message box
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgCont = {'\fontsize{10}There is not {\bf{enough memory}} to compute the {\bf{materialization matrix}} for LoOP!, due to the incorrect input neighborhood parameter {\bf{minPts}}.'; '\fontsize{10}Thus, please {\color{red}\bf{stop}} the execution and run it again with the true optimal parameters, or even just free up some memory space! ;-)'; '\newline{\color{red}\bf{Note:}} Consider the {\bf{BlckSzLim}} value in GUI as a parameter too!'};
            h = msgbox(msgCont,'Memory Error','error',CreateStruct);
            uiwait(h)
            keyboard % Stop the operation from here, because of lack of memory
        end
        %-----------------------------%
        
        matlMatX = [];
        matlMatXind = [];
        cntLim = ceil(n/BLK_SZ_LIM);
        
        for c1 = 1:cntLim
            %     c1
            if c1 ~= cntLim
                Y1 = H.labDS.X((c1-1)*BLK_SZ_LIM+1:c1*BLK_SZ_LIM,:);
            else
                Y1 = H.labDS.X((c1-1)*BLK_SZ_LIM+1:end,:);
            end
            chnkSz = size(Y1,1);
            matlMatY = [];
            matlMatYind = [];
            
            for c2 = 1:cntLim
                %         c2
                if c2 ~= cntLim
                    indVec = (c2-1)*BLK_SZ_LIM+1:c2*BLK_SZ_LIM;
                else
                    indVec = (c2-1)*BLK_SZ_LIM+1:n;
                end
                Y2 = H.labDS.X(indVec,:);
                [sDistY, sDistYind] = sort([pdist2(Y1,Y2) matlMatY],2);
                
                if c1 ~= c2
                    kDist = sDistY(:,kVal);
                    kNNlastInd = find(any(sDistY<=kDist),1,'last');
                    matlMatY = sDistY(:,1:kNNlastInd);
                    
                    truInd = [repmat(indVec,chnkSz,1) matlMatYind];
                    sDistYlinInd = (sDistYind-1).*chnkSz+[1:chnkSz]';
                    matlMatYind = truInd(sDistYlinInd);
                    matlMatYind = matlMatYind(:,1:kNNlastInd);
                else
                    kDist = sDistY(:,kVal+1);
                    kNNlastInd = find(any(sDistY<=kDist),1,'last');
                    matlMatY = sDistY(:,2:kNNlastInd);
                    
                    truInd = [repmat(indVec,chnkSz,1) matlMatYind];
                    sDistYlinInd = (sDistYind-1).*chnkSz+[1:chnkSz]';
                    matlMatYind = truInd(sDistYlinInd);
                    matlMatYind = matlMatYind(:,2:kNNlastInd);
                    
                end
                
                H.progLevl_statText.String = [num2str((c1-1)*cntLim+c2) '/' num2str(cntLim^2)]; pause(.001);
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

