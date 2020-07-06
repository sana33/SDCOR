
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ceit.aut.ac.ir/~sann_cv/
% June 2020

function varargout = MAIN(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAIN_OpeningFcn, ...
                   'gui_OutputFcn',  @MAIN_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function MAIN_OpeningFcn(hO, eventdata, H, varargin)

H.output = hO;

warning off;
PCMact(hO,eventdata,H);

% Update H structure
guidata(hO, H);

function varargout = MAIN_OutputFcn(hO,eventdata,H) 

varargout{1} = H.output;

function load_pushBtn_Callback(hO,eventdata,H)

H.dispOn = get(H.dispOn_checkBox,'Value');

[FileName,PathName] = uigetfile('*.mat', 'Select the dataset along with outlier labels, all as a single MAT-file','..\datasets\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    labDS = importdata([PathName FileName]);
    H.DS = labDS.X;
    H.labFin = labDS.y;
    clearvars labDS
    
    H.OLno = sum(H.labFin==1);
    [H.n,H.p] = size(H.DS);
    
    [~,fileName,~] = fileparts(FileName);
    H.dsName = [fileName '_(' num2str(H.n) 'by' num2str(H.p) ')'];
    H.dsName_statText.String = H.dsName;
    
    if H.dispOn
        if H.p>2
            [coeff_DS_PCA,DS_PCA,~] = pca(H.DS);
            H.coef_PCA = coeff_DS_PCA(:,1:2);
            H.DS_PCA = DS_PCA(:,1:2);
            H.xLim = [min(H.DS_PCA(:,1))-1, max(H.DS_PCA(:,1))+1];
            H.yLim = [min(H.DS_PCA(:,2))-1, max(H.DS_PCA(:,2))+1];
            plotOptional(H,{},'loadDS');
        else
            H.coef_PCA = [];
            H.DS_PCA = [];
            H.xLim = [min(H.DS(:,1))-1, max(H.DS(:,1))+1];
            H.yLim = [min(H.DS(:,2))-1, max(H.DS(:,2))+1];
            plotOptional(H,{},'loadDS');
        end
    else
        H.coef_PCA = [];
        H.DS_PCA = [];
        H.xLim = [min(H.DS(:,1))-1, max(H.DS(:,1))+1];
        H.yLim = [min(H.DS(:,2))-1, max(H.DS(:,2))+1];
    end
    
    msgbox('File was loaded successfully!','Success');
end

guidata(hO,H);

function start_pushBtn_Callback(hO,eventdata,H)

global BLK_SZ_LIM

hOact(hO,eventdata,H,1);
H.dispOn = get(H.dispOn_checkBox,'Value');
H.auxiFig_checkBox.Value = 0;

%------- error handing -------%
if ~isfield(H,'DS') || isempty(H.DS)
    errordlg('Dataset file not found! Please load the input data first!','File Error');
    return
end
if H.dispOn && H.p>2 && isempty(H.coef_PCA)
    [coeff_DS_PCA,DS_PCA,~] = pca(H.DS);
    H.coef_PCA = coeff_DS_PCA(:,1:2);
    H.DS_PCA = DS_PCA(:,1:2);
    H.xLim = [min(H.DS_PCA(:,1))-1, max(H.DS_PCA(:,1))+1];
    H.yLim = [min(H.DS_PCA(:,2))-1, max(H.DS_PCA(:,2))+1];
end
%-----------------------------%

H.dsName_statText.String = H.dsName;
H.chunkSz = str2double(get(H.chunkSz_editText,'String'));
H.PCvarRat = str2double(get(H.PCvarRat_editText,'String'))/100;
H.alphaMemb = str2double(get(H.alphaMemb_editText,'String'));
H.betaPrun = str2double(get(H.betaPrun_editText,'String'));
H.sampRate = str2double(get(H.sampRate_editText,'String'))/100;
H.scorDSszCoef = str2double(get(H.scorDSszCoef_editText,'String'));
BLK_SZ_LIM = str2double(get(H.blckSzlim_editText,'String'));

H.minPtsIntv = str2num(get(H.minPtsIntv_editText,'String'));
H.kStepLngth = str2double(get(H.kStepLngth_editText,'String'));
H.LoOPlambda = str2double(get(H.LoOPlambda_editText,'String'));

H.PCM = get(get(H.PCM_radioBtnGroup,'SelectedObject'),'tag');
H.dimCoef = str2double(get(H.dimCoef_editText,'String'));
H.PSO_particleNo = str2double(get(H.particleNo_editText,'String'));
H.PSO_maxIter = str2double(get(H.maxIter_editText,'String'));
H.PSO_W = str2double(get(H.W_editText,'String'));
H.PSO_C1 = str2double(get(H.C1_editText,'String'));
H.PSO_C2 = str2double(get(H.C2_editText,'String'));
H.PSO_alpha = str2double(get(H.alpha_editText,'String'));
H.manuEps = str2double(get(H.manuEps_editText,'String'));
H.manuMnPt = str2double(get(H.manuMnPt_editText,'String'));
H.epsCoeff = str2double(get(H.epsCoef_editText,'String'));
H.MinPtsCoeff = str2double(get(H.MinPtsCoef_editText,'String'));

if get(H.LOF_checkBox,'Value')
    H.apprType = 'LOF';
    [H.lofVals,H.lofKmat,H.finalAUC,H.tElapsed] = LOF(H);
    set(H.finalAUCbyScores_statText,'String',num2str(H.finalAUC,'%0.3f'));
    set(H.runTime_statText,'String',num2str(H.tElapsed,'%0.3f'));
    msgbox('Process was conducted successfully!','Success');
    
    hOact(hO,eventdata,H,0);
elseif get(H.LoOP_checkBox,'Value')
    H.apprType = 'LoOP';
    [H.LoOPvals,H.finalAUC,H.tElapsed] = LoOP(H);
    set(H.finalAUCbyScores_statText,'String',num2str(H.finalAUC,'%0.3f'));
    set(H.runTime_statText,'String',num2str(H.tElapsed,'%0.3f'));
    msgbox('Process was conducted successfully!','Success');
    
    hOact(hO,eventdata,H,0);
else
    H.apprType = 'SDCOR';
    SDCOR(hO,H);
    H = guidata(hO);
    set(H.runTime_statText,'String',num2str(H.tElapsed,'%0.3f'));
    msgbox('Process was conducted successfully!','Success');
    
    hOact(hO,eventdata,H,0);
    PCMact(hO,eventdata,H);
end

guidata(hO,H);

function topNols_editText_Callback(hO,eventdata,H)

function topNols_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function alphaMemb_editText_Callback(hO,eventdata,H)

function alphaMemb_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function PCvarRat_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function sampRate_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function clearAxes_pushBtn_Callback(hO,eventdata,H)

cla(H.axes1); legend(H.axes1,'off');
cla(H.axes2); legend(H.axes2,'off');
cla(H.axes3); legend(H.axes3,'off');

H.dsName_statText.String = '';
H.finalAUCbyScores_statText.String = '';
H.runTime_statText.String = '';
H.origK_val_statText.String = '';

guidata(hO,H);

function sampRate_editText_Callback(hO,eventdata,H)

function PCvarRat_editText_Callback(hO,eventdata,H)

function PCvarRat_editText_KeyPressFcn(hO,eventdata,H)

function sampRate_editText_KeyPressFcn(hO,eventdata,H)

function finalAUCbyScores_statText_CreateFcn(hO,eventdata,H)

function plotTopNols_pushBtn_Callback(hO,eventdata,H)

if isfield(H,'mahalScores')
    topNolNo = str2double(get(H.topNols_editText,'String'));
    [~,scrSortInd] = sort(H.mahalScores,'descend');
    H.topNol = zeros(H.n,1);
    H.topNol(scrSortInd(1:topNolNo)) = 1;
    
    plotOptional(H,{},'topNol');
else
    msgbox('Sorry! Nothing has been run to display!','Failure','error');
end

function chunkSz_editText_Callback(hO,eventdata,H)

function chunkSz_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function plotLabDS_pushBtn_Callback(hO,eventdata,H)

plotOptional(H,{},'loadDS');

function sampDS_pushBtn_Callback(hO,eventdata,H)

plotOptional(H,{},'sampDS');

function retSetinRed_pushBtn_Callback(hO,eventdata,H)

plotOptional(H,{},'retSetDS');

function scrDS_pushBtn_Callback(hO,eventdata,H)

H.scorDSszCoef = str2double(get(H.scorDSszCoef_editText,'String'));
plotOptional(H,{},'scorDS');

function dimCoef_editText_Callback(hO,eventdata,H)

function dimCoef_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function particleNo_editText_Callback(hO,eventdata,H)

function particleNo_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function maxIter_editText_Callback(hO,eventdata,H)

function maxIter_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function W_editText_Callback(hO,eventdata,H)

function W_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function C1_editText_Callback(hO,eventdata,H)

function C1_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function C2_editText_Callback(hO,eventdata,H)

function C2_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function alpha_editText_Callback(hO,eventdata,H)

function alpha_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function finalMeans_pushBtn_Callback(hO,eventdata,H)

plotOptional(H,{},'finalMeans');

guidata(hO,H);

function regenDS_pushBtn_Callback(hO,eventdata,H)

plotOptional(H,{},'regenDS');

guidata(hO,H);

function auxiFig_checkBox_Callback(hO,eventdata,H)

function epsCoef_editText_Callback(hO,eventdata,H)

function epsCoef_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function MinPtsCoef_editText_Callback(hO,eventdata,H)

function MinPtsCoef_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function manuEps_editText_Callback(hO,eventdata,H)

function manuEps_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function manuMnPt_editText_Callback(hO,eventdata,H)

function manuMnPt_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function scorDSszCoef_editText_Callback(hO,eventdata,H)

function scorDSszCoef_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function saveWork_pushBtn_Callback(hO,eventdata,H)

if isfield(H,'finalAUC')
    [~,fileName,~] = fileparts(H.dsName);
    Hsave = saveWork(hO,eventdata,H);
    
    switch H.apprType
        case 'SDCOR'
            resType = 'SDCOR(visRAM)_';
        case 'LOF'
            resType = 'LOF(visRAM)_';
        case 'LoOP'
            resType = 'LoOP(visRAM)_';
    end
    
    uisave({'Hsave'},['..\results\',resType,'result_$',fileName,'$_AUC=',num2str(H.finalAUC,'%0.3f'),'_Time=',num2str(H.tElapsed,'%0.3f')]);
    msgbox('File was saved successfully!','Success');
else
    msgbox('Sorry! There is not any clear run inf. to be saved!','Failure','error');
end

function loadWork_pushBtn_Callback(hO,eventdata,H)

[FileName,PathName] = uigetfile('*.mat', 'Select the saved workspace to be loaded','..\results\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    Hsave = importdata([PathName FileName]);
    loadWork(hO, eventdata, H, Hsave);
    msgbox('File was loaded successfully!','Success');
end

function [Hsave] = saveWork(hO,eventdata,H)

global BLK_SZ_LIM

switch H.apprType
    case 'SDCOR'
        Hsave = struct('dsName',H.dsName,'labFin',H.labFin,'OLno',H.OLno,'DS_PCA',H.DS_PCA,'coef_PCA',H.coef_PCA,...
            'n',H.n,'p',H.p,'xLim',H.xLim,'yLim',H.yLim,'apprType',H.apprType,'dispOn',H.dispOn,...
            'chunkSz',H.chunkSz,'PCvarRat',str2double(get(H.PCvarRat_editText,'String')),'alphaMemb',H.alphaMemb,'betaPrun',H.betaPrun,...
            'sampRate',str2double(get(H.sampRate_editText,'String')),'topNols',str2double(H.topNols_editText.String),...
            'scorDSszCoef',H.scorDSszCoef,'BLK_SZ_LIM',BLK_SZ_LIM,'PCM',H.PCM,'dimCoef',H.dimCoef,'PSO_particleNo',H.PSO_particleNo,...
            'PSO_maxIter',H.PSO_maxIter,'PSO_W',H.PSO_W,'PSO_C1',H.PSO_C1,'PSO_C2',H.PSO_C2,'PSO_alpha',H.PSO_alpha,...
            'manuEps',H.manuEps,'manuMnPt',H.manuMnPt,'epsCoeff',H.epsCoeff,'MinPtsCoeff',H.MinPtsCoeff,'sampIdx',H.sampIdx,'sampData',H.sampData,...
            'sampData_PCA',H.sampData_PCA,'paramSampDS',H.paramSampDS,'idxSamp',H.idxSamp,...
            'paramCostArrSamp',{H.paramCostArrSamp},'epsilonFin',H.epsilonFin,'MinPtsFin',H.MinPtsFin,...
            'accResArr',H.accResArr,'clusts',{H.clusts},'retIdx',H.retIdx,'means',H.means,'means_PCA',H.means_PCA,'idxMeans',H.idxMeans,...
            'finalClusts',{H.finalClusts},'origK',H.origK,'meansMeans',H.meansMeans,'meansMeans_PCA',H.meansMeans_PCA,...
            'regenDS',H.regenDS,'regenDS_PCA',H.regenDS_PCA,'idxRegenDS',H.idxRegenDS,'origKvec',H.origKvec,...
            'mahalScores',H.mahalScores,'idxFin',H.idxFin,'finalAUC',H.finalAUC,'tElapsed',H.tElapsed);
    case 'LOF'
        Hsave = struct('apprType',H.apprType,'dsName',H.dsName,'labFin',H.labFin,'BLK_SZ_LIM',BLK_SZ_LIM,'minPtsIntv',H.minPtsIntv,...
            'kStepLngth',H.kStepLngth,'lofVals',H.lofVals,'lofKmat',H.lofKmat,'finalAUC',H.finalAUC,'tElapsed',H.tElapsed);
    case 'LoOP'
        Hsave = struct('apprType',H.apprType,'dsName',H.dsName,'labFin',H.labFin,'BLK_SZ_LIM',BLK_SZ_LIM,'minPtsIntv',H.minPtsIntv,...
            'kStepLngth',H.kStepLngth,'LoOPlambda',H.LoOPlambda,'LoOPvals',H.LoOPvals,'finalAUC',H.finalAUC,'tElapsed',H.tElapsed);
end

qOptions.Interpreter = 'tex';
qOptions.Default = 'No';
qstring = 'Would you like to save the input dataset along with the other objects?';
saveChoice = questdlg(qstring,'Save Option','Yes','No',qOptions);
switch saveChoice
    case 'Yes'
        Hsave.DS = H.DS;
        Hsave = orderfields(Hsave);
    case 'No'
        Hsave.DS = [];
        Hsave = orderfields(Hsave);
    case ''
        Hsave.DS = [];
        Hsave = orderfields(Hsave);
end

function [] = loadWork(hO, eventdata, H, Hsave)

global BLK_SZ_LIM

switch Hsave.apprType
    case 'SDCOR'
        H.apprType = Hsave.apprType;
        hOact(hO, eventdata, H, 0);
        H.dsName_statText.String = Hsave.dsName;
        H.dsName = Hsave.dsName;
        H.labFin = Hsave.labFin;
        H.OLno = Hsave.OLno;
        H.DS = Hsave.DS;
        if isempty(Hsave.DS)
            h = msgbox('The dataset file should be loaded separately!','Warning!','warn');
            uiwait(h)
        end
        H.DS_PCA = Hsave.DS_PCA;
        H.coef_PCA = Hsave.coef_PCA;
        H.n = Hsave.n;
        H.p = Hsave.p;
        H.xLim = Hsave.xLim;
        H.yLim = Hsave.yLim;
        H.dispOn = Hsave.dispOn;
        if Hsave.dispOn
            H.dispOn_checkBox.Value = 1;
            plotOptional(H,{},'loadDS');
        else
            H.dispOn_checkBox.Value = 0;
        end
        
        H.chunkSz_editText.String = num2str(Hsave.chunkSz);
        H.PCvarRat_editText.String = num2str(Hsave.PCvarRat);
        H.alphaMemb_editText.String = num2str(Hsave.alphaMemb);
        H.betaPrun_editText.String = num2str(Hsave.betaPrun);
        H.sampRate_editText.String = num2str(Hsave.sampRate);
        H.topNols_editText.String = num2str(Hsave.topNols);
        H.scorDSszCoef_editText.String = num2str(Hsave.scorDSszCoef);
        BLK_SZ_LIM = Hsave.BLK_SZ_LIM;
        H.blckSzlim_editText.String = num2str(BLK_SZ_LIM);
        
        switch Hsave.PCM
            case 'PSO_pcm_radioBtn'
                H.PSO_pcm_radioBtn.Value = 1;
                
                H.PSO_finCond = 1;
                plotOptional(H,{Hsave.paramCostArrSamp{1},Hsave.paramSampDS},'PSOcost');
                
            case 'manu_pcm_radioBtn'
                H.manu_pcm_radioBtn.Value = 1;
                
                cla(H.axes3);
        end
        PCMact(hO,eventdata,H);
        
        H.accResArr = Hsave.accResArr;
        H.accFinCond = 1;
        plotOptional(H,{Hsave.accResArr},'accPerChunk');
        
        H.dimCoef_editText.String = num2str(Hsave.dimCoef);
        H.particleNo_editText.String = num2str(Hsave.PSO_particleNo);
        H.maxIter_editText.String = num2str(Hsave.PSO_maxIter);
        H.W_editText.String = num2str(Hsave.PSO_W);
        H.C1_editText.String = num2str(Hsave.PSO_C1);
        H.C2_editText.String = num2str(Hsave.PSO_C2);
        H.alpha_editText.String = num2str(Hsave.PSO_alpha);
        H.manuEps_editText.String = num2str(Hsave.manuEps);
        H.manuMnPt_editText.String = num2str(Hsave.manuMnPt);
        H.epsCoef_editText.String = num2str(Hsave.epsCoeff);
        H.MinPtsCoef_editText.String = num2str(Hsave.MinPtsCoeff);
        
        H.sampIdx = Hsave.sampIdx;
        H.sampData = Hsave.sampData;
        H.sampData_PCA = Hsave.sampData_PCA;
        H.paramCostArrSamp = Hsave.paramCostArrSamp;
        H.paramSampDS = Hsave.paramSampDS;
        H.idxSamp = Hsave.idxSamp;
        H.epsilonFin = Hsave.epsilonFin;
        H.MinPtsFin = Hsave.MinPtsFin;
        
        H.clusts = Hsave.clusts;
        H.retIdx = Hsave.retIdx;
        H.means = Hsave.means;
        H.means_PCA = Hsave.means_PCA;
        H.idxMeans = Hsave.idxMeans;
        
        H.finalClusts = Hsave.finalClusts;
        H.origK = Hsave.origK;
        H.meansMeans = Hsave.meansMeans;
        H.meansMeans_PCA = Hsave.meansMeans_PCA;
        H.regenDS = Hsave.regenDS;
        H.regenDS_PCA = Hsave.regenDS_PCA;
        H.idxRegenDS = Hsave.idxRegenDS;
        H.origKvec = Hsave.origKvec;
        
        H.mahalScores = Hsave.mahalScores;
        H.idxFin = Hsave.idxFin;
        H.finalAUC = Hsave.finalAUC;
        H.tElapsed = Hsave.tElapsed;
        H.finalAUCbyScores_statText.String = num2str(H.finalAUC,'%0.3f');
        H.runTime_statText.String = num2str(Hsave.tElapsed,'%0.3f');
        H.origK_val_statText.String = num2str(Hsave.origK);
        
    case 'LOF'
        H.apprType = Hsave.apprType;
        H.LOF_checkBox.Value = 1;
        hOact(hO, eventdata, H, 0);
        
        H.dsName = Hsave.dsName;
        H.dsName_statText.String = Hsave.dsName;
        H.labFin = Hsave.labFin;
        H.DS = Hsave.DS;
        if isempty(Hsave.DS)
            h = msgbox('The dataset file should be loaded separately!','Warning!','warn');
            uiwait(h)
        end
        BLK_SZ_LIM = Hsave.BLK_SZ_LIM;
        H.blckSzlim_editText.String = num2str(BLK_SZ_LIM);
        
        H.minPtsIntv = Hsave.minPtsIntv;
        H.kStepLngth = Hsave.kStepLngth;
        H.lofVals = Hsave.lofVals;
        H.lofKmat = Hsave.lofKmat;
        H.finalAUC = Hsave.finalAUC;
        H.tElapsed = Hsave.tElapsed;
        H.finalAUCbyScores_statText.String = num2str(H.finalAUC,'%0.3f');
        H.runTime_statText.String = num2str(Hsave.tElapsed);

    case 'LoOP'
        H.apprType = Hsave.apprType;
        H.LoOP_checkBox.Value = 1;
        hOact(hO, eventdata, H, 0);
        
        H.dsName = Hsave.dsName;
        H.dsName_statText.String = Hsave.dsName;
        H.labFin = Hsave.labFin;
        H.DS = Hsave.DS;
        if isempty(Hsave.DS)
            h = msgbox('The dataset file should be loaded separately!','Warning!','warn');
            uiwait(h)
        end
        BLK_SZ_LIM = Hsave.BLK_SZ_LIM;
        H.blckSzlim_editText.String = num2str(BLK_SZ_LIM);
        
        H.minPtsIntv = Hsave.minPtsIntv;
        H.kStepLngth = Hsave.kStepLngth;
        H.LoOPlambda = Hsave.LoOPlambda;
        H.LoOPvals = Hsave.LoOPvals;
        H.finalAUC = Hsave.finalAUC;
        H.tElapsed = Hsave.tElapsed;
        H.finalAUCbyScores_statText.String = num2str(H.finalAUC,'%0.3f');
        H.runTime_statText.String = num2str(Hsave.tElapsed);

end

guidata(hO,H);

function PCMact(hO,eventdata,H)

switch H.PCM_radioBtnGroup.SelectedObject.Tag
    case 'PSO_pcm_radioBtn'
        H.makeManu_pushBtn.Enable =  'on';
        
        H.dimCoef_editText.Enable = 'on';
        H.particleNo_editText.Enable = 'on';
        H.maxIter_editText.Enable = 'on';
        H.W_editText.Enable = 'on';
        H.C1_editText.Enable = 'on';
        H.C2_editText.Enable = 'on';
        H.alpha_editText.Enable = 'on';
        H.manuEps_editText.Enable = 'off';
        H.manuMnPt_editText.Enable = 'off';
        H.epsCoef_editText.Enable = 'on';
        H.MinPtsCoef_editText.Enable = 'on';
        
    case 'manu_pcm_radioBtn'
        H.makeManu_pushBtn.Enable =  'off';
        
        H.dimCoef_editText.Enable = 'off';
        H.particleNo_editText.Enable = 'off';
        H.maxIter_editText.Enable = 'off';
        H.W_editText.Enable = 'off';
        H.C1_editText.Enable = 'off';
        H.C2_editText.Enable = 'off';
        H.alpha_editText.Enable = 'off';
        H.manuEps_editText.Enable = 'on';
        H.manuMnPt_editText.Enable = 'on';
        H.epsCoef_editText.Enable = 'on';
        H.MinPtsCoef_editText.Enable = 'on';
        
end

function hOact(hO, eventdata, H, startCond)

if startCond
    SDCOR_InitParam_Act(hO,eventdata,H,0);
    clearAxes_pushBtn_Callback(hO,eventdata,H);
    
    H.LOF_checkBox.Enable = 'off';
    H.LoOP_checkBox.Enable = 'off';
    LOF_LoOP_InitParam_Act(hO,eventdata,H,0);
    LOF_InitParam_Act(hO,eventdata,H,0);
    LoOP_InitParam_Act(hO,eventdata,H,0);
        
    dispPlot_Act(hO,eventdata,H,0);
    
    mainButns_Act(hO,eventdata,H,0);
    
    DBSCANparamCM_Act(hO,eventdata,H,0);

else
    if strcmp(H.apprType,'SDCOR')
        SDCOR_InitParam_Act(hO,eventdata,H,1);
        
        H.LOF_checkBox.Enable = 'on';
        H.LoOP_checkBox.Enable = 'on';
        H.LOF_checkBox.Value = 0;
        H.LoOP_checkBox.Value = 0;
        LOF_LoOP_InitParam_Act(hO,eventdata,H,1);
        LOF_InitParam_Act(hO,eventdata,H,1);
        LoOP_InitParam_Act(hO,eventdata,H,1);
        
        dispPlot_Act(hO,eventdata,H,1);
        
        DBSCANparamCM_Act(hO,eventdata,H,1);
        
    elseif strcmp(H.apprType,'LOF')
        H.LOF_checkBox.Enable = 'on';
        LOF_checkBox_Callback(hO,eventdata,H);
        
    elseif strcmp(H.apprType,'LoOP')
        H.LoOP_checkBox.Enable = 'on';
        LoOP_checkBox_Callback(hO,eventdata,H);
        
    end
    
    mainButns_Act(hO,eventdata,H,1);
end

function SDCOR_InitParam_Act(hO, eventdata, H, actCond)

if ~actCond
    H.chunkSz_editText.Enable = 'off';
    H.PCvarRat_editText.Enable = 'off';
    H.alphaMemb_editText.Enable = 'off';
    H.betaPrun_editText.Enable = 'off';
    H.sampRate_editText.Enable = 'off';
    H.topNols_editText.Enable = 'off';
    H.scorDSszCoef_editText.Enable = 'off';
    H.blckSzlim_editText.Enable = 'off';
    
    H.dispOn_checkBox.Enable = 'off';
else
    H.chunkSz_editText.Enable = 'on';
    H.PCvarRat_editText.Enable = 'on';
    H.alphaMemb_editText.Enable = 'on';
    H.betaPrun_editText.Enable = 'on';
    H.sampRate_editText.Enable = 'on';
    H.topNols_editText.Enable = 'on';
    H.scorDSszCoef_editText.Enable = 'on';
    H.blckSzlim_editText.Enable = 'on';
    
    H.dispOn_checkBox.Enable = 'on';
end

function dispPlot_Act(hO, eventdata, H, actCond)

if ~actCond
    H.auxiFig_checkBox.Enable = 'off';
    H.plotLabDS_pushBtn.Enable = 'off';
    H.sampDS_pushBtn.Enable = 'off';
    H.retSetinRed_pushBtn.Enable = 'off';
    H.finalMeans_pushBtn.Enable = 'off';
    H.regenDS_pushBtn.Enable = 'off';
    H.scrDS_pushBtn.Enable = 'off';
    H.plotTopNols_pushBtn.Enable = 'off';
else
    H.auxiFig_checkBox.Enable = 'on';
    H.plotLabDS_pushBtn.Enable = 'on';
    H.sampDS_pushBtn.Enable = 'on';
    H.retSetinRed_pushBtn.Enable = 'on';
    H.finalMeans_pushBtn.Enable = 'on';
    H.regenDS_pushBtn.Enable = 'on';
    H.scrDS_pushBtn.Enable = 'on';
    H.plotTopNols_pushBtn.Enable = 'on';
end

function mainButns_Act(hO, eventdata, H, actCond)

if ~actCond
    H.load_pushBtn.Enable = 'off';
    H.start_pushBtn.Enable = 'off';
    H.clearAxes_pushBtn.Enable = 'off';
    H.saveWork_pushBtn.Enable = 'off';
    H.loadWork_pushBtn.Enable = 'off';
else
    H.load_pushBtn.Enable = 'on';
    H.start_pushBtn.Enable = 'on';
    H.clearAxes_pushBtn.Enable = 'on';
    H.saveWork_pushBtn.Enable = 'on';
    H.loadWork_pushBtn.Enable = 'on';
end

function DBSCANparamCM_Act(hO, eventdata, H, actCond)

if ~actCond
    H.PSO_pcm_radioBtn.Enable = 'off';
    H.manu_pcm_radioBtn.Enable = 'off';
    H.makeManu_pushBtn.Enable =  'off';
    
    H.dimCoef_editText.Enable = 'off';
    H.particleNo_editText.Enable = 'off';
    H.maxIter_editText.Enable = 'off';
    H.W_editText.Enable = 'off';
    H.C1_editText.Enable = 'off';
    H.C2_editText.Enable = 'off';
    H.alpha_editText.Enable = 'off';
    H.manuEps_editText.Enable = 'off';
    H.manuMnPt_editText.Enable = 'off';
    H.epsCoef_editText.Enable = 'off';
    H.MinPtsCoef_editText.Enable = 'off';
else
    H.PSO_pcm_radioBtn.Enable = 'on';
    H.manu_pcm_radioBtn.Enable = 'on';
    H.makeManu_pushBtn.Enable =  'on';
    
    H.dimCoef_editText.Enable = 'on';
    H.particleNo_editText.Enable = 'on';
    H.maxIter_editText.Enable = 'on';
    H.W_editText.Enable = 'on';
    H.C1_editText.Enable = 'on';
    H.C2_editText.Enable = 'on';
    H.alpha_editText.Enable = 'on';
    H.manuEps_editText.Enable = 'on';
    H.manuMnPt_editText.Enable = 'on';
    H.epsCoef_editText.Enable = 'on';
    H.MinPtsCoef_editText.Enable = 'on';
end
   
function PSO_pcm_radioBtn_Callback(hO,eventdata,H)

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
msgCont = ['\fontsize{10}For excluding the execution time of PSO, after finding the optimal parameters, ' ...
    'press the {\bf{Make Manu}} button to set them as manual, and then run the algorithm again.'];
h = msgbox(msgCont,'Memory Error','error',CreateStruct);
uiwait(h)

PCMact(hO,eventdata,H);

function manu_pcm_radioBtn_Callback(hO,eventdata,H)

PCMact(hO,eventdata,H);

function makeManu_pushBtn_Callback(hO,eventdata,H)

if isfield(H,'paramSampDS')
    H.manuEps_editText.String = num2str(H.paramSampDS(1));
    H.manuMnPt_editText.String = num2str(H.paramSampDS(2));
else
    msgbox('Sorry! Nothing has been run to set for!','Failure','error');
end

function resetBtns_pushBtn_Callback(hO,eventdata,H)

H.apprType = 'SDCOR';
hOact(hO, eventdata, H, 0);
PCMact(hO,eventdata,H);

function LOF_checkBox_Callback(hO,eventdata,H)

if get(H.LOF_checkBox,'Value')
    H.LoOP_checkBox.Enable = 'off';
    H.LoOP_checkBox.Value = 0;
    LOF_LoOP_InitParam_Act(hO,eventdata,H,1);
    LOF_InitParam_Act(hO,eventdata,H,1);
    LoOP_InitParam_Act(hO,eventdata,H,0);
    
    SDCOR_InitParam_Act(hO,eventdata,H,0);
    H.blckSzlim_editText.Enable = 'on';
    
    dispPlot_Act(hO,eventdata,H,0);
    H.dispOn_checkBox.Value = 0;
    H.auxiFig_checkBox.Value = 0;
    
    DBSCANparamCM_Act(hO,eventdata,H,0);
else
    H.LoOP_checkBox.Enable = 'on';
    LOF_LoOP_InitParam_Act(hO,eventdata,H,1);
    LOF_InitParam_Act(hO,eventdata,H,1);
    LoOP_InitParam_Act(hO,eventdata,H,1);
    
    SDCOR_InitParam_Act(hO,eventdata,H,1);
    
    dispPlot_Act(hO,eventdata,H,1);
    
    DBSCANparamCM_Act(hO,eventdata,H,1);
    PCMact(hO,eventdata,H);
end

function LoOP_checkBox_Callback(hO,eventdata,H)

if get(H.LoOP_checkBox,'Value')
    H.LOF_checkBox.Enable = 'off';
    H.LOF_checkBox.Value = 0;
    LOF_LoOP_InitParam_Act(hO,eventdata,H,1);
    LOF_InitParam_Act(hO,eventdata,H,0);
    LoOP_InitParam_Act(hO,eventdata,H,1);
    
    H.minPtsIntv_Temp = H.minPtsIntv_editText.String;
    H.minPtsIntv_editText.String = num2str(mean(str2num(H.minPtsIntv_editText.String)));
    
    SDCOR_InitParam_Act(hO,eventdata,H,0);
    H.blckSzlim_editText.Enable = 'on';
    
    dispPlot_Act(hO,eventdata,H,0);
    H.dispOn_checkBox.Value = 0;
    H.auxiFig_checkBox.Value = 0;
    
    DBSCANparamCM_Act(hO,eventdata,H,0);
else
    H.LOF_checkBox.Enable = 'on';
    LOF_LoOP_InitParam_Act(hO,eventdata,H,1);
    LOF_InitParam_Act(hO,eventdata,H,1);
    LoOP_InitParam_Act(hO,eventdata,H,1);
    
    H.minPtsIntv_editText.String = H.minPtsIntv_Temp;
    
    SDCOR_InitParam_Act(hO,eventdata,H,1);
    
    dispPlot_Act(hO,eventdata,H,1);
    
    DBSCANparamCM_Act(hO,eventdata,H,1);
    PCMact(hO,eventdata,H);
end

guidata(hO,H);

function LOF_LoOP_InitParam_Act(hO, eventdata, H, actCond)

if ~actCond
    H.minPtsIntv_editText.Enable = 'off';
else
    H.minPtsIntv_editText.Enable = 'on';
end

function LOF_InitParam_Act(hO, eventdata, H, actCond)

if ~actCond
    H.kStepLngth_editText.Enable = 'off';
else
    H.kStepLngth_editText.Enable = 'on';
end
    
function LoOP_InitParam_Act(hO, eventdata, H, actCond)

if ~actCond
    H.LoOPlambda_editText.Enable = 'off';
else
    H.LoOPlambda_editText.Enable = 'on';
end

function blckSzlim_editText_Callback(hO,eventdata,H)

function blckSzlim_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function dispOn_checkBox_Callback(hO,eventdata,H)

function betaPrun_editText_Callback(hO,eventdata,H)

function betaPrun_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function minPtsIntv_editText_Callback(hO,eventdata,H)

function minPtsIntv_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function kStepLngth_editText_Callback(hO,eventdata,H)

function kStepLngth_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function LoOPlambda_editText_Callback(hO,eventdata,H)

function LoOPlambda_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end
