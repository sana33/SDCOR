
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ce.aut.ac.ir/~sann_cv/
% June 2020

function plotOptional(handles,data,option)

    if isfield(handles,'DS') && isfield(handles,'dispOn')
        switch option
            case 'loadDS'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                contour_levels = 20;
                if handles.p>2 && handles.dispOn && ~isempty(handles.DS_PCA)
                    gscatter(handles.DS_PCA(:,1),handles.DS_PCA(:,2),handles.labFin,'br'); 
                    xlim(handles.xLim); ylim(handles.yLim); grid on; hold on;
                    plotDensityContour(handles.DS_PCA, contour_levels); hold off
                    legend('inliers','outliers','PDF','location','best');
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.DS)
                    gscatter(handles.DS(:,1),handles.DS(:,2),handles.labFin,'br'); 
                    xlim(handles.xLim); ylim(handles.yLim); grid on; hold on;
                    plotDensityContour(handles.DS, contour_levels); hold off
                    legend('inliers','outliers','PDF','location','best');
                elseif ~handles.dispOn || isempty(handles.DS_PCA) || isempty(handles.DS)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'sampDS'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.sampData_PCA)
                    gscatter(handles.sampData_PCA(:,1),handles.sampData_PCA(:,2),handles.idxSamp);
                    xlim(handles.xLim); ylim(handles.yLim); grid on; pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.sampData)
                    gscatter(handles.sampData(:,1),handles.sampData(:,2),handles.idxSamp);
                    xlim(handles.xLim); ylim(handles.yLim); grid on; pause(1);
                elseif ~handles.dispOn || isempty(handles.sampData_PCA) || isempty(handles.sampData)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'retSetDS'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.DS_PCA) && ~isempty(handles.means_PCA)
                    if ~handles.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                    gscatter(handles.DS_PCA(:,1),handles.DS_PCA(:,2),handles.labFin,'br');
                    grid on; hold on;
                    plot(handles.means_PCA(:,1),handles.means_PCA(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                        'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                    plot(handles.DS_PCA(handles.retIdx,1),handles.DS_PCA(handles.retIdx,2),'Marker','o','MarkerSize',3, ...
                        'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                    hold off;
                    legend('inliers','outliers',sprintf('means=%d',size(handles.means,1)),...
                        sprintf('retainSet=%d',numel(handles.retIdx)),'location','best');
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.DS) && ~isempty(handles.means)
                    gscatter(handles.DS(:,1),handles.DS(:,2),handles.labFin,'br');
                    grid on; hold on;
                    plot(handles.means(:,1),handles.means(:,2),'Marker','s','MarkerSize',7,'MarkerFaceColor','m', ...
                        'MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                    plot(handles.DS(handles.retIdx,1),handles.DS(handles.retIdx,2),'Marker','o','MarkerSize',3, ...
                        'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',3,'LineStyle','none');
                    hold off;
                    legend('inliers','outliers',sprintf('means=%d',size(handles.means,1)),...
                        sprintf('retainSet=%d',numel(handles.retIdx)),'location','best');
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif ~handles.dispOn || isempty(handles.DS_PCA) || isempty(handles.DS) || isempty(handles.means_PCA) || isempty(handles.means)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'finalMeans'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.means_PCA) && ~isempty(handles.meansMeans_PCA) && ~isempty(handles.origK)
                    if ~handles.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                    [handles.xLim,handles.yLim] = xyLimCreat(handles.means_PCA,handles.xLim,handles.yLim);
                    cMap = hsv(handles.origK);
                    gscatter(handles.means_PCA(:,1),handles.means_PCA(:,2),handles.idxMeans,cMap);
                    grid on; hold on;
                    hm1 = gscatter(handles.meansMeans_PCA(:,1),handles.meansMeans_PCA(:,2),handles.origKvec,cMap,'^',8);
                    for p1 = 1:handles.origK
                        hm1(p1).MarkerFaceColor = 'k';
                    end
                    hold off;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.means) && ~isempty(handles.meansMeans) && ~isempty(handles.origK)
                    cMap = hsv(handles.origK);
                    gscatter(handles.means(:,1),handles.means(:,2),handles.idxMeans,cMap);
                    grid on; hold on;
                    hm1 = gscatter(handles.meansMeans(:,1),handles.meansMeans(:,2),handles.origKvec,cMap,'^',8);
                    for p1 = 1:handles.origK
                        hm1(p1).MarkerFaceColor = 'k';
                    end
                    hold off;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif ~handles.dispOn || isempty(handles.means_PCA) || isempty(handles.means) || isempty(handles.meansMeans_PCA) || isempty(handles.meansMeans) || ...
                        isempty(handles.origK)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'regenDS'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.regenDS_PCA) && ~isempty(handles.meansMeans_PCA) && ~isempty(handles.origK)
                    if ~handles.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                    [handles.xLim,handles.yLim] = xyLimCreat(handles.regenDS_PCA(:,1:2),handles.xLim,handles.yLim);
                    cMap = hsv(handles.origK);
                    gscatter(handles.regenDS_PCA(:,1),handles.regenDS_PCA(:,2),handles.idxRegenDS,cMap);
                    grid on; hold on;
                    hm1 = gscatter(handles.meansMeans_PCA(:,1),handles.meansMeans_PCA(:,2),handles.origKvec,cMap,'^',8);
                    for p1 = 1:handles.origK
                        hm1(p1).MarkerFaceColor = 'k';
                    end
                    hold off;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.regenDS) && ~isempty(handles.meansMeans) && ~isempty(handles.origK)
                    cMap = hsv(handles.origK);
                    gscatter(handles.regenDS(:,1),handles.regenDS(:,2),handles.idxRegenDS,cMap);
                    grid on; hold on;
                    hm1 = gscatter(handles.meansMeans(:,1),handles.meansMeans(:,2),handles.origKvec,cMap,'^',8);
                    for p1 = 1:handles.origK
                        hm1(p1).MarkerFaceColor = 'k';
                    end
                    hold off;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif ~handles.dispOn || isempty(handles.regenDS_PCA) || isempty(handles.regenDS) || isempty(handles.meansMeans_PCA) || ...
                        isempty(handles.meansMeans) || isempty(handles.origK)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'scorDS'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.DS_PCA) && ~isempty(handles.origK)
                    if ~handles.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                    cMap = hsv(handles.origK);
                    hm1 = gscatter(handles.DS_PCA(:,1),handles.DS_PCA(:,2),handles.idxFin,cMap(handles.idxFin,:),'.',handles.mahalScores.*handles.scorDSszCoef,'on');
                    %             hm1 = gscatter(handles.DS_PCA(:,1),handles.DS_PCA(:,2),handles.idxFin,cMap(handles.idxFin,:),'.',ones(handles.n,1).*handles.scorDSszCoef,'on');
                    %             for p1 = 1:handles.origK
                    %                 hm1(p1).MarkerFaceColor = cMap(p1,:);
                    %                 hm1(p1).MarkerEdgeColor = cMap(p1,:);
                    %             end
                    grid on;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.DS) && ~isempty(handles.origK)
                    cMap = hsv(handles.origK);
                    hm1 = gscatter(handles.DS(:,1),handles.DS(:,2),handles.idxFin,cMap(handles.idxFin,:),'.',handles.mahalScores.*handles.scorDSszCoef,'on');
                    %             hm1 = gscatter(handles.DS(:,1),handles.DS(:,2),handles.idxFin,cMap(handles.idxFin,:),'.',ones(handles.n,1).*handles.scorDSszCoef,'on');
                    %             for p1 = 1:handles.origK
                    %                 hm1(p1).MarkerFaceColor = cMap(p1,:);
                    %                 hm1(p1).MarkerEdgeColor = cMap(p1,:);
                    %             end
                    grid on;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif ~handles.dispOn || isempty(handles.DS_PCA) || isempty(handles.DS) || isempty(handles.origK)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'topNol'
                if ~get(handles.auxiFig_checkBox,'Value'); axes(handles.axes1); else; figure; end
                if handles.p>2 && handles.dispOn && ~isempty(handles.DS_PCA) && ~isempty(handles.topNol)
                    if ~handles.dispOn; msgbox('Sorry! Nothing to display for p>2!','Failure','error'); end
                    gscatter(handles.DS_PCA(:,1),handles.DS_PCA(:,2),handles.topNol,'br','.');
                    legend('inliers','outliers'); grid on;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif handles.p==2 && handles.dispOn && ~isempty(handles.DS) && ~isempty(handles.topNol)
                    gscatter(handles.DS(:,1),handles.DS(:,2),handles.topNol,'br','.');
                    legend('inliers','outliers'); grid on;
                    xlim(handles.xLim); ylim(handles.yLim);
                    pause(1);
                elseif ~handles.dispOn || isempty(handles.DS_PCA) || isempty(handles.DS) || isempty(handles.topNol)
                    uiwait(msgbox(['Sorry! Required items (maybe including the query dataset in 2D case) have not been saved to display. '...
                        'Please preform a fresh run with the "DispPlot" option checked, for detailed illustrations.'],'Failure','error','modal'));
                end
                
            case 'accPerChunk'
                axes(handles.axes2);
                xAxis = 1:size(data{1},2);
                plot(xAxis,data{1}(1,:),'-r',xAxis,data{1}(2,:),'-k',xAxis,data{1}(3,:),'-g',xAxis,data{1}(4,:),'-b',xAxis,data{1}(4,:),'-m');
                pause(.001);
                if handles.accFinCond
                    legend(sprintf('ROC per chunk\nFinal Value=%0.3f',full(data{1}(1,end))), ...
                        sprintf('PR per chunk\nFinal Value=%0.3f',full(data{1}(2,end))), ...
                        sprintf('Precision per chunk\nFinal Value=%0.3f',full(data{1}(3,end))), ...
                        sprintf('Recall per chunk\nFinal Value=%0.3f',full(data{1}(4,end))), ...
                        sprintf('F1-Score per chunk\nFinal Value=%0.3f',full(data{1}(5,end))),'location','best');
                    grid on;
                end
                
            case 'PSOcost'
                axes(handles.axes3);
                plot(data{1},'-r');
                pause(.001);
                if handles.PSO_finCond
                    legend(sprintf('PSO costArr for SampDS\nEps=%0.3f, MinPts=%d',full(data{2}(1)),full(data{2}(2))),'location','best');
                    grid on;
                end
                
        end
    else
        uiwait(msgbox('Sorry! Nothing has been run to display.','Failure','error','modal'));
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

function plotDensityContour(data, contour_levels)
    % Check if the input data is not empty
    if isempty(data)
        error('Input data cannot be empty');
    end
    
    % Create a grid over the data region
    x_min = min(data(:, 1)) - 1;  % Add some padding
    x_max = max(data(:, 1)) + 1;
    y_min = min(data(:, 2)) - 1;
    y_max = max(data(:, 2)) + 1;
    
    % Define grid resolution
    [X, Y] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));
    
    % Estimate the probability density function using ksdensity
    [density, ~] = ksdensity(data, [X(:), Y(:)]);
    
    % Reshape the density values to match the grid shape
    Z = reshape(density, size(X));
    
    % Plot contour lines of different levels for the PDF
    contour(X, Y, Z, contour_levels, 'DisplayName', 'PDF');  % contour_levels levels of contours
    colorbar; 
    
end
