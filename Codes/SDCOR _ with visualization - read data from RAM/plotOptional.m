function plotOptional(H,data,option)

if isfield(H,'p') && isfield(H,'dispOn')
    switch option
        case 'loadDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.labFin,'gb');
            elseif H.p<=2 && H.dispOn
                gscatter(H.DS(:,1),H.DS(:,2),H.labFin,'gb');
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            grid on;
            xlim(H.xLim); ylim(H.yLim);
            
        case 'sampDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                gscatter(H.sampData_PCA(:,1),H.sampData_PCA(:,2),H.idxSamp);
            elseif H.p<=2 && H.dispOn
                gscatter(H.sampData(:,1),H.sampData(:,2),H.idxSamp);
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            grid on;
            xlim(H.xLim); ylim(H.yLim);
            pause(1);
            
        case 'retSetDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                if ~H.dispOn; h = msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.labFin,'gb');
                grid on; hold on;
                plot(H.means_PCA(:,1),H.means_PCA(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                    'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                plot(H.DS_PCA(H.retIdx,1),H.DS_PCA(H.retIdx,2),'Marker','o','MarkerSize',3, ...
                    'MarkerFaceColor','r','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                hold off;
                legend('inliers','outliers',sprintf('means=%d',size(H.means,1)),...
                    sprintf('retainSet=%d',numel(H.retIdx)),'location','best');
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p<=2 && H.dispOn
                gscatter(H.DS(:,1),H.DS(:,2),H.labFin,'gb');
                grid on; hold on;
                plot(H.means(:,1),H.means(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                    'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                plot(H.DS(H.retIdx,1),H.DS(H.retIdx,2),'Marker','o','MarkerSize',3, ...
                    'MarkerFaceColor','r','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                hold off;
                legend('inliers','outliers',sprintf('means=%d',size(H.means,1)),...
                    sprintf('retainSet=%d',numel(H.retIdx)),'location','best');
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            
        case 'finalMeans'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                if ~H.dispOn; h = msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
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
            elseif H.p<=2 && H.dispOn
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
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            
        case 'regenDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                if ~H.dispOn; h = msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
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
            elseif H.p<=2 && H.dispOn
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
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            
        case 'scorDS'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                if ~H.dispOn; h = msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
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
            elseif H.p<=2 && H.dispOn
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
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            
        case 'topNol'
            if ~get(H.auxiFig_checkBox,'Value'); axes(H.axes1); else; figure; end
            if H.p>2 && H.dispOn
                if ~H.dispOn; h = msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                gscatter(H.DS_PCA(:,1),H.DS_PCA(:,2),H.topNol,'br','.',[5 10]);
                grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif H.p<=2 && H.dispOn
                gscatter(H.DS(:,1),H.DS(:,2),H.topNol,'br','.',[5 10]);
                grid on;
                xlim(H.xLim); ylim(H.yLim);
                pause(1);
            elseif ~H.dispOn
                h = msgbox('Sorry! Nothing has been saved to display!','Failure','error');
            end
            
        case 'accPerChunk'
            axes(H.axes2);
            xAxis = 1:size(data{1},2);
            plot(xAxis,data{1}(1,:),'-r',xAxis,data{1}(2,:),'-k',xAxis,data{1}(3,:),'-g',xAxis,data{1}(4,:),'-b');
            pause(.001);
            if H.accFinCond
                legend(sprintf('AUC per chunk\nFinal Value=%0.2f',data{1}(1,end)), ...
                    sprintf('Precision per chunk\nFinal Value=%0.2f',data{1}(2,end)), ...
                    sprintf('Recall per chunk\nFinal Value=%0.2f',data{1}(3,end)), ...
                    sprintf('F1-Score per chunk\nFinal Value=%0.2f',data{1}(4,end)),'location','best');
                grid on;
            end
            
        case 'PSOcost'
            axes(H.axes3);
            plot(data{1},'-r');
            pause(.001);
            if H.PSO_finCond
                legend(sprintf('PSO costArr for SampDS\neps=%0.2f, MinPts=%d', ...
                    data{2}(1),data{2}(2)), ...
                    'location','best');
                grid on;
            end
            
    end
else
    h = msgbox('Sorry! Nothing has been run to display!','Failure','error');
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

