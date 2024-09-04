
% Author: Sayyed-Ahmad Naghavi-Nozad, M.Sc., Artificial Intelligence
% AmirKabir University of Technology, Department of Computer Engineering
% Email Address: sa_na33@aut.ac.ir, ahmad.naghavi.aut@gmail.com
% Website: https://ce.aut.ac.ir/~sann_cv/
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

function MAIN_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

warning off;
handles.startCond = 0;
PCMact(hObject,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

function varargout = MAIN_OutputFcn(hObject,eventdata,handles) 

varargout{1} = handles.output;

function load_pushBtn_Callback(hObject,eventdata,handles)

clearAxes_pushBtn_Callback(hObject,eventdata,handles); % clearing workspace before loading new data
handles.dispOn = get(handles.dispOn_checkBox,'Value');

[FileName,PathName] = uigetfile('*.mat', 'Select the dataset along with outlier labels, all as a single MAT-file','..\datasets\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    Hclear(hObject,eventdata,handles); handles = guidata(hObject);
    labDS = importdata([PathName FileName]);
    handles.DS = labDS.X;
    handles.labFin = labDS.y;
    clearvars labDS
    
    handles.OLno = sum(handles.labFin==1);
    [handles.n,handles.p] = size(handles.DS);
    
    [~,handles.dsName,~] = fileparts(FileName);
    handles.dsName_statText.String = handles.dsName;
    handles.chunkSz_editText.String = ceil(.1*handles.n);
    handles.topNols_editText.String = floor(.03*handles.n);
    
    if handles.dispOn
        if handles.p>2
            [coeff_DS_PCA,DS_PCA,~] = pca(handles.DS);
            handles.coef_PCA = coeff_DS_PCA(:,1:2);
            handles.DS_PCA = DS_PCA(:,1:2);
            handles.xLim = [min(handles.DS_PCA(:,1))-1, max(handles.DS_PCA(:,1))+1];
            handles.yLim = [min(handles.DS_PCA(:,2))-1, max(handles.DS_PCA(:,2))+1];
            plotOptional(handles,{},'loadDS');
        else
            handles.coef_PCA = [];
            handles.DS_PCA = [];
            handles.xLim = [min(handles.DS(:,1))-1, max(handles.DS(:,1))+1];
            handles.yLim = [min(handles.DS(:,2))-1, max(handles.DS(:,2))+1];
            plotOptional(handles,{},'loadDS');
        end
    else
        handles.coef_PCA = [];
        handles.DS_PCA = [];
        handles.xLim = [min(handles.DS(:,1))-1, max(handles.DS(:,1))+1];
        handles.yLim = [min(handles.DS(:,2))-1, max(handles.DS(:,2))+1];
    end
    
    handles.manu_pcm_radioBtn.Value = 1; PCMact(hObject,eventdata,handles);
    msgbox('File was loaded successfully!','Success');
end

guidata(hObject,handles);

function start_pushBtn_Callback(hObject,eventdata,handles)

handles.startCond = 1; hOact(hObject,eventdata,handles);
handles.dispOn = get(handles.dispOn_checkBox,'Value');
handles.auxiFig_checkBox.Value = 0;

%------- error handing -------%
if ~isfield(handles,'DS') || isempty(handles.DS)
    errordlg('Dataset file not found! Please load the input data first!','File Error');
    return
end

if handles.dispOn && handles.p>2 && isempty(handles.coef_PCA)
    [coeff_DS_PCA,DS_PCA,~] = pca(handles.DS);
    handles.coef_PCA = coeff_DS_PCA(:,1:2);
    handles.DS_PCA = DS_PCA(:,1:2);
    handles.xLim = [min(handles.DS_PCA(:,1))-1, max(handles.DS_PCA(:,1))+1];
    handles.yLim = [min(handles.DS_PCA(:,2))-1, max(handles.DS_PCA(:,2))+1];
end
%-----------------------------%

handles.dsName_statText.String = handles.dsName;
handles.chunkSz = str2double(get(handles.chunkSz_editText,'String'));
handles.PCvarRat = str2double(get(handles.PCvarRat_editText,'String'))/100;
handles.alphaMemb = str2double(get(handles.alphaMemb_editText,'String'));
handles.betaPrun = str2double(get(handles.betaPrun_editText,'String'));
handles.sampRate = str2double(get(handles.sampRate_editText,'String'))/100;
handles.scorDSszCoef = str2double(get(handles.scorDSszCoef_editText,'String'));
handles.nonUnifSamp = get(handles.nonUnifSamp_chckbx,'Value');

handles.PCM = get(get(handles.PCM_radioBtnGroup,'SelectedObject'),'tag');
handles.PSO_particleNo = str2double(get(handles.particleNo_editText,'String'));
handles.PSO_maxIter = str2double(get(handles.maxIter_editText,'String'));
handles.PSO_W = str2double(get(handles.W_editText,'String'));
handles.PSO_C1 = str2double(get(handles.C1_editText,'String'));
handles.PSO_C2 = str2double(get(handles.C2_editText,'String'));
handles.PSO_alpha = str2double(get(handles.alpha_editText,'String'));
handles.manuEps = str2double(get(handles.manuEps_editText,'String'));
if isnan(handles.manuEps); handles.manuEps = 0; end
handles.manuMnPt = str2double(get(handles.manuMnPt_editText,'String'));
if isnan(handles.manuMnPt); handles.manuMnPt = floor(log(handles.n)); handles.manuMnPt_editText.String = num2str(handles.manuMnPt); end
handles.epsCoeff = str2double(get(handles.epsCoef_editText,'String'));

SDCOR(hObject,handles);
handles = guidata(hObject);
set(handles.finalROC_statText,'String',num2str(handles.finalROC,'%0.3f'));
set(handles.finalPR_statText,'String',num2str(handles.finalPR,'%0.3f'));
set(handles.runTime_statText,'String',num2str(handles.tElapsed,'%0.3f'));
msgbox('Process was conducted successfully!','Success');

handles.startCond = 0; hOact(hObject,eventdata,handles);
PCMact(hObject,eventdata,handles);

guidata(hObject,handles);

function topNols_editText_Callback(hObject,eventdata,handles)

function topNols_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alphaMemb_editText_Callback(hObject,eventdata,handles)

function alphaMemb_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PCvarRat_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sampRate_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clearAxes_pushBtn_Callback(hObject,eventdata,handles)

cla(handles.axes1); legend(handles.axes1,'off');
cla(handles.axes2); legend(handles.axes2,'off');
cla(handles.axes3); legend(handles.axes3,'off');

handles.dsName_statText.String = '';
handles.finalROC_statText.String = '';
handles.finalPR_statText.String = '';
handles.runTime_statText.String = '';
handles.origK_val_statText.String = '';
if handles.startCond==0; handles.manuMnPt_editText.String = ''; end

guidata(hObject,handles);

function Hclear(hObject,eventdata,handles)

Harr = {'DS','dsName','labFin','OLno','DS_PCA','coef_PCA','n','p','xLim','yLim','chunkSz','PCvarRat','alphaMemb',...
    'betaPrun','sampRate','topNols','scorDSszCoef','PCM','PSO_particleNo','PSO_maxIter','PSO_W',...
    'PSO_C1','PSO_C2','PSO_alpha','manuEps','manuMnPt','epsCoeff','paramSampDS','paramCostArrSamp','origEps','origMnPt','origK',...
    'sampInd','sampData','sampData_PCA','idxSamp','accResArr','clusts','retIdx','means','means_PCA','idxMeans','finalClusts',...
    'meansMeans','meansMeans_PCA','regenDS','regenDS_PCA','idxRegenDS','origKvec','mahalScores','idxFin','finalROC','finalPR','tElapsed'};

for c1 = 1:numel(Harr)
    if isfield(handles,Harr{c1})
        handles = setfield(handles,Harr{c1},[]);
    end
end

guidata(hObject,handles);

function sampRate_editText_Callback(hObject,eventdata,handles)

function PCvarRat_editText_Callback(hObject,eventdata,handles)

function PCvarRat_editText_KeyPressFcn(hObject,eventdata,handles)

function sampRate_editText_KeyPressFcn(hObject,eventdata,handles)

function finalROC_statText_CreateFcn(hObject,eventdata,handles)

function chunkSz_editText_Callback(hObject,eventdata,handles)

function chunkSz_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotLabDS_pushBtn_Callback(hObject,eventdata,handles)

plotOptional(handles,{},'loadDS');

function sampDS_pushBtn_Callback(hObject,eventdata,handles)

plotOptional(handles,{},'sampDS');

function retSetinRed_pushBtn_Callback(hObject,eventdata,handles)

plotOptional(handles,{},'retSetDS');

function finalMeans_pushBtn_Callback(hObject,eventdata,handles)

plotOptional(handles,{},'finalMeans');

guidata(hObject,handles);

function regenDS_pushBtn_Callback(hObject,eventdata,handles)

plotOptional(handles,{},'regenDS');

guidata(hObject,handles);

function scrDS_pushBtn_Callback(hObject,eventdata,handles)

handles.scorDSszCoef = str2double(get(handles.scorDSszCoef_editText,'String'));
plotOptional(handles,{},'scorDS');

function plotTopNols_pushBtn_Callback(hObject,eventdata,handles)

if isfield(handles,'mahalScores') && ~isempty(handles.mahalScores)
    topNolNo = str2double(get(handles.topNols_editText,'String')); topNolNo(isnan(topNolNo)) = 0;
    [~,scrSortInd] = sort(handles.mahalScores,'descend');
    handles.topNol = zeros(handles.n,1);
    handles.topNol(scrSortInd(1:topNolNo)) = 1;
else
    handles.topNol = [];
end
plotOptional(handles,{},'topNol');

function particleNo_editText_Callback(hObject,eventdata,handles)

function particleNo_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxIter_editText_Callback(hObject,eventdata,handles)

function maxIter_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function W_editText_Callback(hObject,eventdata,handles)

function W_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function C1_editText_Callback(hObject,eventdata,handles)

function C1_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function C2_editText_Callback(hObject,eventdata,handles)

function C2_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alpha_editText_Callback(hObject,eventdata,handles)

function alpha_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function auxiFig_checkBox_Callback(hObject,eventdata,handles)

function epsCoef_editText_Callback(hObject,eventdata,handles)

function epsCoef_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function epsCoef_editText_ButtonDownFcn(hObject, eventdata, handles)

function manuEps_editText_Callback(hObject,eventdata,handles)

function manuEps_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function manuMnPt_editText_Callback(hObject,eventdata,handles)

function manuMnPt_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scorDSszCoef_editText_Callback(hObject,eventdata,handles)

function scorDSszCoef_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function saveWork_pushBtn_Callback(hObject,eventdata,handles)

if isfield(handles,'finalROC') && ~isempty(handles.finalROC)
    Hsave = saveWork(hObject,eventdata,handles);
    
    resType = 'SDCOR(visRAM)_';
    uisave({'Hsave'},['..\results\',resType,'result_$',handles.dsName,'$_ROC=',num2str(handles.finalROC,'%0.3f'),'_PR=',num2str(handles.finalPR,'%0.3f'),...
        '_Time=',num2str(handles.tElapsed,'%0.3f'),'.mat']);
    msgbox('File was saved successfully!','Success');
else
    msgbox('Sorry! There is not any clear run info to be saved!','Failure','error');
end

function loadWork_pushBtn_Callback(hObject,eventdata,handles)

[FileName,PathName] = uigetfile('*.mat', 'Select the saved workspace to be loaded','..\results\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    Hsave = importdata([PathName FileName]);
    loadWork(hObject, eventdata, handles, Hsave);
    msgbox('File was loaded successfully!','Success');
end

function [Hsave] = saveWork(hObject,eventdata,handles)

Hsave = struct('dsName',handles.dsName,'labFin',handles.labFin,'OLno',handles.OLno,'DS_PCA',handles.DS_PCA,'coef_PCA',handles.coef_PCA,'n',handles.n,'p',handles.p,...
    'xLim',handles.xLim,'yLim',handles.yLim,'dispOn',handles.dispOn,'chunkSz',handles.chunkSz,'PCvarRat',handles.PCvarRat*100,'alphaMemb',handles.alphaMemb,...
    'betaPrun',handles.betaPrun,'sampRate',handles.sampRate*100,'topNols',str2double(handles.topNols_editText.String),'scorDSszCoef',handles.scorDSszCoef,...
    'nonUnifSamp',handles.nonUnifSamp,'PCM',handles.PCM,'PSO_particleNo',handles.PSO_particleNo,'PSO_maxIter',handles.PSO_maxIter,'PSO_W',handles.PSO_W,...
    'PSO_C1',handles.PSO_C1,'PSO_C2',handles.PSO_C2,'PSO_alpha',handles.PSO_alpha,'manuEps',handles.manuEps,'manuMnPt',handles.manuMnPt,'epsCoeff',handles.epsCoeff,...
    'paramSampDS',handles.paramSampDS,'paramCostArrSamp',{handles.paramCostArrSamp},'origEps',handles.origEps,'origMnPt',handles.origMnPt,'origK',handles.origK,...
    'sampInd',handles.sampInd,'sampData',handles.sampData,'sampData_PCA',handles.sampData_PCA,'idxSamp',handles.idxSamp,'accResArr',handles.accResArr,...
    'clusts',{handles.clusts},'retIdx',handles.retIdx,'means',handles.means,'means_PCA',handles.means_PCA,'idxMeans',handles.idxMeans,'finalClusts',{handles.finalClusts},...
    'meansMeans',handles.meansMeans,'meansMeans_PCA',handles.meansMeans_PCA,'regenDS',handles.regenDS,'regenDS_PCA',handles.regenDS_PCA,...
    'idxRegenDS',handles.idxRegenDS,'origKvec',handles.origKvec,'mahalScores',handles.mahalScores,'idxFin',handles.idxFin,'finalROC',handles.finalROC,...
    'finalPR',handles.finalPR,'tElapsed',handles.tElapsed);

qOptions.Interpreter = 'tex';
qOptions.Default = 'No';
qstring = 'Would you like to save the input dataset along with the other objects?';
saveChoice = questdlg(qstring,'Save Option','Yes','No',qOptions);
switch saveChoice
    case 'Yes'
        Hsave.DS = handles.DS;
        Hsave = orderfields(Hsave);
    case 'No'
        Hsave.DS = [];
        Hsave = orderfields(Hsave);
    case ''
        Hsave.DS = [];
        Hsave = orderfields(Hsave);
end

function [] = loadWork(hObject, eventdata, handles, Hsave)

clearAxes_pushBtn_Callback(hObject,eventdata,handles); % clearing workspace before loading the saved result
Hclear(hObject,eventdata,handles); handles = guidata(hObject);

handles.startCond = 0; hOact(hObject,eventdata,handles);
handles.dsName_statText.String = Hsave.dsName; handles.dsName = Hsave.dsName;
handles.dsName = Hsave.dsName;
handles.labFin = Hsave.labFin;
handles.OLno = Hsave.OLno;
handles.DS = Hsave.DS;
handles.DS_PCA = Hsave.DS_PCA;
handles.coef_PCA = Hsave.coef_PCA;
handles.n = Hsave.n;
handles.p = Hsave.p;
handles.xLim = Hsave.xLim;
handles.yLim = Hsave.yLim;
handles.dispOn = Hsave.dispOn;
if Hsave.dispOn && ~isempty(Hsave.DS)
    handles.dispOn_checkBox.Value = 1;
    plotOptional(handles,{},'loadDS');
elseif ~Hsave.dispOn
    handles.dispOn_checkBox.Value = 0;
end

handles.chunkSz_editText.String = num2str(Hsave.chunkSz); handles.chunkSz = Hsave.chunkSz;
handles.PCvarRat_editText.String = num2str(Hsave.PCvarRat); handles.PCvarRat = Hsave.PCvarRat/100;
handles.alphaMemb_editText.String = num2str(Hsave.alphaMemb); handles.alphaMemb = Hsave.alphaMemb;
handles.betaPrun_editText.String = num2str(Hsave.betaPrun); handles.betaPrun = Hsave.betaPrun;
handles.sampRate_editText.String = num2str(Hsave.sampRate); handles.sampRate = Hsave.sampRate/100;
handles.topNols_editText.String = num2str(Hsave.topNols); handles.topNols = Hsave.topNols;
handles.scorDSszCoef_editText.String = num2str(Hsave.scorDSszCoef); handles.scorDSszCoef = Hsave.scorDSszCoef;
handles.nonUnifSamp_chckbx.Value = Hsave.nonUnifSamp; handles.nonUnifSamp = Hsave.nonUnifSamp;

handles.PCM = Hsave.PCM;
switch Hsave.PCM
    case 'PSO_pcm_radioBtn'
        handles.PSO_pcm_radioBtn.Value = 1;
        handles.PSO_finCond = 1;
        plotOptional(handles,{Hsave.paramCostArrSamp{1},Hsave.paramSampDS},'PSOcost');
        
    case 'manu_pcm_radioBtn'
        handles.manu_pcm_radioBtn.Value = 1;
        cla(handles.axes3);
        
end
PCMact(hObject,eventdata,handles);

handles.accResArr = Hsave.accResArr;
handles.accFinCond = 1;
plotOptional(handles,{Hsave.accResArr},'accPerChunk');

handles.particleNo_editText.String = num2str(Hsave.PSO_particleNo); handles.PSO_particleNo = Hsave.PSO_particleNo;
handles.maxIter_editText.String = num2str(Hsave.PSO_maxIter); handles.PSO_maxIter = Hsave.PSO_maxIter;
handles.W_editText.String = num2str(Hsave.PSO_W); handles.PSO_W = Hsave.PSO_W;
handles.C1_editText.String = num2str(Hsave.PSO_C1); handles.PSO_C1 = Hsave.PSO_C1;
handles.C2_editText.String = num2str(Hsave.PSO_C2); handles.PSO_C2 = Hsave.PSO_C2;
handles.alpha_editText.String = num2str(Hsave.PSO_alpha); handles.PSO_alpha = Hsave.PSO_alpha;
handles.manuEps_editText.String = num2str(Hsave.manuEps); handles.manuEps = Hsave.manuEps;
handles.manuMnPt_editText.String = num2str(Hsave.manuMnPt); handles.manuMnPt = Hsave.manuMnPt;
handles.epsCoef_editText.String = num2str(Hsave.epsCoeff); handles.epsCoeff = Hsave.epsCoeff;
handles.origK_val_statText.String = num2str(Hsave.origK); handles.origK = Hsave.origK;

handles.sampInd = Hsave.sampInd;
handles.sampData = Hsave.sampData;
handles.sampData_PCA = Hsave.sampData_PCA;
handles.paramCostArrSamp = Hsave.paramCostArrSamp;
handles.paramSampDS = Hsave.paramSampDS;
handles.idxSamp = Hsave.idxSamp;
handles.origEps = Hsave.origEps;
handles.origMnPt = Hsave.origMnPt;

handles.clusts = Hsave.clusts;
handles.retIdx = Hsave.retIdx;
handles.means = Hsave.means;
handles.means_PCA = Hsave.means_PCA;
handles.idxMeans = Hsave.idxMeans;

handles.finalClusts = Hsave.finalClusts;
handles.origK = Hsave.origK;
handles.meansMeans = Hsave.meansMeans;
handles.meansMeans_PCA = Hsave.meansMeans_PCA;
handles.regenDS = Hsave.regenDS;
handles.regenDS_PCA = Hsave.regenDS_PCA;
handles.idxRegenDS = Hsave.idxRegenDS;
handles.origKvec = Hsave.origKvec;

handles.mahalScores = Hsave.mahalScores;
handles.idxFin = Hsave.idxFin;
handles.finalROC_statText.String = num2str(Hsave.finalROC,'%0.3f'); handles.finalROC = Hsave.finalROC;
handles.finalPR_statText.String = num2str(Hsave.finalPR,'%0.3f'); handles.finalPR = Hsave.finalPR;
handles.runTime_statText.String = num2str(Hsave.tElapsed,'%0.3f'); handles.tElapsed = Hsave.tElapsed;

guidata(hObject,handles);

function PCMact(hObject,eventdata,handles)

switch handles.PCM_radioBtnGroup.SelectedObject.Tag
    case 'Kgraph_pcm_radioBtn'
        handles.makeManu_pushBtn.Enable =  'off';
        
        handles.particleNo_editText.Enable = 'off';
        handles.maxIter_editText.Enable = 'off';
        handles.W_editText.Enable = 'off';
        handles.C1_editText.Enable = 'off';
        handles.C2_editText.Enable = 'off';
        handles.alpha_editText.Enable = 'off';
        handles.manuEps_editText.Enable = 'off';
        handles.manuMnPt_editText.Enable = 'on';
        handles.epsCoef_editText.Enable = 'off';
        
    case 'PSO_pcm_radioBtn'
        handles.makeManu_pushBtn.Enable =  'on';
        
        handles.particleNo_editText.Enable = 'on';
        handles.maxIter_editText.Enable = 'on';
        handles.W_editText.Enable = 'on';
        handles.C1_editText.Enable = 'on';
        handles.C2_editText.Enable = 'on';
        handles.alpha_editText.Enable = 'on';
        handles.manuEps_editText.Enable = 'off';
        handles.manuMnPt_editText.Enable = 'on';
        handles.epsCoef_editText.Enable = 'on';
        
    case 'manu_pcm_radioBtn'
        handles.makeManu_pushBtn.Enable =  'off';
        
        handles.particleNo_editText.Enable = 'off';
        handles.maxIter_editText.Enable = 'off';
        handles.W_editText.Enable = 'off';
        handles.C1_editText.Enable = 'off';
        handles.C2_editText.Enable = 'off';
        handles.alpha_editText.Enable = 'off';
        handles.manuEps_editText.Enable = 'on';
        handles.manuMnPt_editText.Enable = 'on';
        handles.epsCoef_editText.Enable = 'on';
        
end

function hOact(hObject, eventdata, handles)

if handles.startCond
    SDCOR_InitParam_Act(hObject,eventdata,handles,0);
    clearAxes_pushBtn_Callback(hObject,eventdata,handles);
    dispPlot_Act(hObject,eventdata,handles,0);
    mainButns_Act(hObject,eventdata,handles,0);
    DBSCANparamCM_Act(hObject,eventdata,handles,0);

else
    SDCOR_InitParam_Act(hObject,eventdata,handles,1);
    dispPlot_Act(hObject,eventdata,handles,1);
    DBSCANparamCM_Act(hObject,eventdata,handles,1);
    mainButns_Act(hObject,eventdata,handles,1);
    
end

function SDCOR_InitParam_Act(hObject, eventdata, handles, actCond)

if ~actCond
    handles.chunkSz_editText.Enable = 'off';
    handles.PCvarRat_editText.Enable = 'off';
    handles.alphaMemb_editText.Enable = 'off';
    handles.betaPrun_editText.Enable = 'off';
    handles.sampRate_editText.Enable = 'off';
    handles.topNols_editText.Enable = 'off';
    handles.scorDSszCoef_editText.Enable = 'off';
    handles.blckSzlim_editText.Enable = 'off';
    
    handles.dispOn_checkBox.Enable = 'off';
else
    handles.chunkSz_editText.Enable = 'on';
    handles.PCvarRat_editText.Enable = 'on';
    handles.alphaMemb_editText.Enable = 'on';
    handles.betaPrun_editText.Enable = 'on';
    handles.sampRate_editText.Enable = 'on';
    handles.topNols_editText.Enable = 'on';
    handles.scorDSszCoef_editText.Enable = 'on';
    handles.blckSzlim_editText.Enable = 'on';
    
    handles.dispOn_checkBox.Enable = 'on';
end

function dispPlot_Act(hObject, eventdata, handles, actCond)

if ~actCond
    handles.auxiFig_checkBox.Enable = 'off';
    handles.plotLabDS_pushBtn.Enable = 'off';
    handles.sampDS_pushBtn.Enable = 'off';
    handles.retSetinRed_pushBtn.Enable = 'off';
    handles.finalMeans_pushBtn.Enable = 'off';
    handles.regenDS_pushBtn.Enable = 'off';
    handles.scrDS_pushBtn.Enable = 'off';
    handles.plotTopNols_pushBtn.Enable = 'off';
else
    handles.auxiFig_checkBox.Enable = 'on';
    handles.plotLabDS_pushBtn.Enable = 'on';
    handles.sampDS_pushBtn.Enable = 'on';
    handles.retSetinRed_pushBtn.Enable = 'on';
    handles.finalMeans_pushBtn.Enable = 'on';
    handles.regenDS_pushBtn.Enable = 'on';
    handles.scrDS_pushBtn.Enable = 'on';
    handles.plotTopNols_pushBtn.Enable = 'on';
end

function mainButns_Act(hObject, eventdata, handles, actCond)

if ~actCond
    handles.load_pushBtn.Enable = 'off';
    handles.start_pushBtn.Enable = 'off';
    handles.clearAxes_pushBtn.Enable = 'off';
    handles.saveWork_pushBtn.Enable = 'off';
    handles.loadWork_pushBtn.Enable = 'off';
else
    handles.load_pushBtn.Enable = 'on';
    handles.start_pushBtn.Enable = 'on';
    handles.clearAxes_pushBtn.Enable = 'on';
    handles.saveWork_pushBtn.Enable = 'on';
    handles.loadWork_pushBtn.Enable = 'on';
end

function DBSCANparamCM_Act(hObject, eventdata, handles, actCond)

if ~actCond
    handles.Kgraph_pcm_radioBtn.Enable = 'off';
    handles.PSO_pcm_radioBtn.Enable = 'off';
    handles.manu_pcm_radioBtn.Enable = 'off';
    handles.makeManu_pushBtn.Enable =  'off';
    
    handles.particleNo_editText.Enable = 'off';
    handles.maxIter_editText.Enable = 'off';
    handles.W_editText.Enable = 'off';
    handles.C1_editText.Enable = 'off';
    handles.C2_editText.Enable = 'off';
    handles.alpha_editText.Enable = 'off';
    handles.manuEps_editText.Enable = 'off';
    handles.manuMnPt_editText.Enable = 'off';
    handles.epsCoef_editText.Enable = 'off';
else
    handles.Kgraph_pcm_radioBtn.Enable = 'on';
    handles.PSO_pcm_radioBtn.Enable = 'on';
    handles.manu_pcm_radioBtn.Enable = 'on';
    handles.makeManu_pushBtn.Enable =  'on';
    
    handles.particleNo_editText.Enable = 'on';
    handles.maxIter_editText.Enable = 'on';
    handles.W_editText.Enable = 'on';
    handles.C1_editText.Enable = 'on';
    handles.C2_editText.Enable = 'on';
    handles.alpha_editText.Enable = 'on';
    handles.manuEps_editText.Enable = 'on';
    handles.manuMnPt_editText.Enable = 'on';
    handles.epsCoef_editText.Enable = 'on';
end

function Kgraph_pcm_radioBtn_Callback(hObject, eventdata, handles)

handles.manuEps_editText.String = '';

if isfield(handles,'p') && ~isempty(handles.p)
    handles.manuMnPt_editText.String = 10*handles.p;
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    msgCont = '\fontsize{10} To avoid outliers, we set the {\it{MinPts}} quantity to {\bf{10{\cdot}p}}. You can change it at your will!';
    uiwait(msgbox(msgCont,'Alert!','Help',CreateStruct));
end

PCMact(hObject,eventdata,handles);
   
function PSO_pcm_radioBtn_Callback(hObject,eventdata,handles)

handles.manuEps_editText.String = '';

if isfield(handles,'p') && ~isempty(handles.p)
    handles.manuMnPt_editText.String = 10*handles.p;
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    msgCont = ['\fontsize{10} PSO searches in predetermined ranges by the user for finding the optimal {\it{Eps}} and {\it{MinPts}} ',...
        'quantities. The {\it{Eps}} range and the lower bound for the {\it{MinPts}} range are defined automatically; for the {\it{MinPts}} ',...
        'upper bound, we set it to {\bf{10{\cdot}p}} for your convenience.\newline\newline {\color{red}\bf{Note:}} For being more prudent ',...
        'in the case of high noisy datasets, you can change it to greater values at your will! But please do NOT go much further as for ',...
        'vary large {\it{MinPts}} values, the following {\it{k}}-dist graph tends to become horizontal and without any distinguishable ',...
        'valleys; therefore, it might not be possible to choose the right value for the {\it{MinPts}} parameter.'];
    uiwait(msgbox(msgCont,'Alert!','Help',CreateStruct));
end

PCMact(hObject,eventdata,handles);

function manu_pcm_radioBtn_Callback(hObject,eventdata,handles)

PCMact(hObject,eventdata,handles);

function makeManu_pushBtn_Callback(hObject,eventdata,handles)

if isfield(handles,'paramSampDS') && ~isempty(handles.paramSampDS)
    handles.manuEps_editText.String = num2str(handles.paramSampDS(1));
    handles.manuMnPt_editText.String = num2str(handles.paramSampDS(2));
else
    msgbox('Sorry! Nothing has been run to set for!','Failure','error');
end

function resetBtns_pushBtn_Callback(hObject,eventdata,handles)

handles.startCond = 0; hOact(hObject, eventdata, handles);
PCMact(hObject,eventdata,handles);

guidata(hObject,handles);

function blckSzlim_editText_Callback(hObject,eventdata,handles)

function blckSzlim_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dispOn_checkBox_Callback(hObject,eventdata,handles)

function betaPrun_editText_Callback(hObject,eventdata,handles)

function betaPrun_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minPtsIntv_editText_Callback(hObject,eventdata,handles)

function minPtsIntv_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kStepLngth_editText_Callback(hObject,eventdata,handles)

function kStepLngth_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nonUnifSamp_chckbx_Callback(hObject, eventdata, handles)
