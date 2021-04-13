
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ce.aut.ac.ir/~sann_cv/
% June 2020

function plotOptional(H,data,option)

if isfield(H,'DS') && all(size(H.DS)==[H.n,H.p]) && isfield(H,'dispOn')
    switch option
        case 'loadDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.DS_PCA)
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.labFin,'br');
                legend('inliers','outliers','location','best');
                xlim(H.xLim); ylim(H.yLim); grid on;
            elseif H.p==2 && H.dispOn && ~isempty(H.DS)
                gscatter(H.DS(:,1),H.DS(:,2),H.labFin,'br');
                legend('inliers','outliers','location','best');
                xlim(H.xLim); ylim(H.yLim); grid on;
            elseif ~H.dispOn || isempty(H.DS_PCA) || isempty(H.DS)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'sampDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.sampData_PCA)
                gscatter(H.sampData_PCA(:,1),H.sampData_PCA(:,2),H.idxSamp);
                xlim(H.xLim); ylim(H.yLim); grid on; pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.sampData)
                gscatter(H.sampData(:,1),H.sampData(:,2),H.idxSamp);
                xlim(H.xLim); ylim(H.yLim); grid on; pause(1);
            elseif ~H.dispOn || isempty(H.sampData_PCA) || isempty(H.sampData)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'retSetDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.DS_PCA)
                if ~H.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.labFin,'br');
                grid on; hold on;
                plot(H.means_PCA(:,1),H.means_PCA(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                    'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                plot(H.DS_PCA(H.retIdx,1),H.DS_PCA(H.retIdx,2),'Marker','o','MarkerSize',3, ...
                    'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                hold off;
                legend('inliers','outliers',sprintf('means=%d',size(H.means,1)),...
                    sprintf('retainSet=%d',numel(H.retIdx)),'location','best');
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.DS)
                gscatter(H.DS(:,1),H.DS(:,2),H.labFin,'br');
                grid on; hold on;
                plot(H.means(:,1),H.means(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                    'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                plot(H.DS(H.retIdx,1),H.DS(H.retIdx,2),'Marker','o','MarkerSize',3, ...
                    'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                hold off;
                legend('inliers','outliers',sprintf('means=%d',size(H.means,1)),...
                    sprintf('retainSet=%d',numel(H.retIdx)),'location','best');
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn || isempty(H.DS_PCA) || isempty(H.DS)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'finalMeans'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.means_PCA)
                if ~H.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                [H.xLim,H.yLim] = xyLimCreat(H.means_PCA,H.xLim,H.yLim);
                cMap = hsv(H.origK);
                gscatter(H.means_PCA(:,1),H.means_PCA(:,2),H.idxMeans,cMap);
                grid on; hold on;
                hm1 = gscatter(H.meansMeans_PCA(:,1),H.meansMeans_PCA(:,2),H.origKvec,cMap,'^',8);
                for p1 = 1:H.origK
                    hm1(p1).MarkerFaceColor = 'k';
                end
                hold off;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.means)
                cMap = hsv(H.origK);
                gscatter(H.means(:,1),H.means(:,2),H.idxMeans,cMap);
                grid on; hold on;
                hm1 = gscatter(H.meansMeans(:,1),H.meansMeans(:,2),H.origKvec,cMap,'^',8);
                for p1 = 1:H.origK
                    hm1(p1).MarkerFaceColor = 'k';
                end
                hold off;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn || isempty(H.means_PCA) || isempty(H.means)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'regenDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.regenDS_PCA)
                if ~H.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                [H.xLim,H.yLim] = xyLimCreat(H.regenDS_PCA(:,1:2),H.xLim,H.yLim);
                cMap = hsv(H.origK);
                gscatter(H.regenDS_PCA(:,1),H.regenDS_PCA(:,2),H.idxRegenDS,cMap);
                grid on; hold on;
                hm1 = gscatter(H.meansMeans_PCA(:,1),H.meansMeans_PCA(:,2),H.origKvec,cMap,'^',8);
                for p1 = 1:H.origK
                    hm1(p1).MarkerFaceColor = 'k';
                end
                hold off;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.regenDS)
                cMap = hsv(H.origK);
                gscatter(H.regenDS(:,1),H.regenDS(:,2),H.idxRegenDS,cMap);
                grid on; hold on;
                hm1 = gscatter(H.meansMeans(:,1),H.meansMeans(:,2),H.origKvec,cMap,'^',8);
                for p1 = 1:H.origK
                    hm1(p1).MarkerFaceColor = 'k';
                end
                hold off;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn || isempty(H.regenDS_PCA) || isempty(H.regenDS)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'scorDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.DS_PCA)
                if ~H.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                cMap = hsv(H.origK);
                hm1 = gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.idxFin,cMap(H.idxFin,:),'.',H.mahalScores.*H.scorDSszCoef,'on');
                %             hm1 = gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.idxFin,cMap(H.idxFin,:),'.',ones(H.n,1).*H.scorDSszCoef,'on');
                %             for p1 = 1:H.origK
                %                 hm1(p1).MarkerFaceColor = cMap(p1,:);
                %                 hm1(p1).MarkerEdgeColor = cMap(p1,:);
                %             end
                grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.DS)
                cMap = hsv(H.origK);
                hm1 = gscatter(H.DS(:,1),H.DS(:,2),H.idxFin,cMap(H.idxFin,:),'.',H.mahalScores.*H.scorDSszCoef,'on');
                %             hm1 = gscatter(H.DS(:,1),H.DS(:,2),H.idxFin,cMap(H.idxFin,:),'.',ones(H.n,1).*H.scorDSszCoef,'on');
                %             for p1 = 1:H.origK
                %                 hm1(p1).MarkerFaceColor = cMap(p1,:);
                %                 hm1(p1).MarkerEdgeColor = cMap(p1,:);
                %             end
                grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn || isempty(H.DS_PCA) || isempty(H.DS)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'topNol'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn && ~isempty(H.DS_PCA)
                if ~H.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.topNol,'br','.');
                legend('inliers','outliers'); grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p==2 && H.dispOn && ~isempty(H.DS)
                gscatter(H.DS(:,1),H.DS(:,2),H.topNol,'br','.');
                legend('inliers','outliers'); grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn || isempty(H.DS_PCA) || isempty(H.DS)
                uiwait(msgbox('Sorry! Required items have not been saved to display!','Failure','error','modal'));
            end
            
        case 'accPerChunk'
            axes(H.axes2);
            xAxis = 1:size(data{1},2);
            plot(xAxis,data{1}(1,:),'-r',xAxis,data{1}(2,:),'-k',xAxis,data{1}(3,:),'-g',xAxis,data{1}(4,:),'-b',xAxis,data{1}(4,:),'-m');
            pause(.001);
            if H.accFinCond
                legend(sprintf('ROC per chunk\nFinal Value=%0.3f',data{1}(1,end)), ...
                    sprintf('PR per chunk\nFinal Value=%0.3f',data{1}(2,end)), ...
                    sprintf('Precision per chunk\nFinal Value=%0.3f',data{1}(3,end)), ...
                    sprintf('Recall per chunk\nFinal Value=%0.3f',data{1}(4,end)), ...
                    sprintf('F1-Score per chunk\nFinal Value=%0.3f',data{1}(5,end)),'location','best');
                grid on;
            end
            
        case 'PSOcost'
            axes(H.axes3);
            plot(data{1},'-r');
            pause(.001);
            if H.PSO_finCond
                legend(sprintf('PSO costArr for SampDS\nEps=%0.3f, MinPts=%d',data{2}(1),data{2}(2)),'location','best');
                grid on;
            end
            
    end
else
    uiwait(msgbox('Sorry! Nothing has been run to display or the query dataset has not been saved for illustration.',...
        'Failure','error','modal'));
end

end

function [xLim,yLim] = xyLimCreat(cents,xLim,yLim)

K = size(cents,1);
if K~=1
    centsMin = min(cents);
    centsMax = max(cents);
else
    centsMin = cents;
    centsMax = cents;
end

xLim = [min([xLim, centsMin(1)])-1, max([xLim, centsMax(1)])+1];
yLim = [min([yLim, centsMin(2)])-1, max([yLim, centsMax(2)])+1];

end

