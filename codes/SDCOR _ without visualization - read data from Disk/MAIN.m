
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
PCMact(hObject,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

function varargout = MAIN_OutputFcn(hObject,eventdata,handles) 

varargout{1} = handles.output;

function load_pushBtn_Callback(hObject,eventdata,handles)

clearAxes_pushBtn_Callback(hObject,eventdata,handles); % clearing workspace before loading new data

[FileName,PathName] = uigetfile('*.mat', 'Select the dataset along with outlier labels, all as a single MAT-file','..\datasets\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    Hclear(hObject,eventdata,handles); handles = guidata(hObject);
    handles.labDS = matfile([PathName FileName]);
    [handles.n,handles.p] = size(handles.labDS,'X');
    
    [~,handles.dsName,~] = fileparts(FileName);
    handles.dsName_statText.String = handles.dsName;
    handles.chunkSz_editText.String = ceil(.1*handles.n);
    
    handles.manu_pcm_radioBtn.Value = 1; PCMact(hObject,eventdata,handles);
    msgbox('File was loaded successfully!','Success');
end

guidata(hObject,handles);

function start_pushBtn_Callback(hObject,eventdata,handles)

%------- error handing -------%
if ~isfield(handles,'labDS') || isempty(handles.labDS)
    errordlg('Dataset file not found! Please load the input data first!','File Error');
    return
end
%-----------------------------%

handles.startCond = 1; hOact(hObject,eventdata,handles);

handles.dsName_statText.String = handles.dsName;
handles.chunkSz = str2double(get(handles.chunkSz_editText,'String'));
handles.PCvarRat = str2double(get(handles.PCvarRat_editText,'String'))/100;
handles.alphaMemb = str2double(get(handles.alphaMemb_editText,'String'));
handles.betaPrun = str2double(get(handles.betaPrun_editText,'String'));
handles.sampRate = str2double(get(handles.sampRate_editText,'String'))/100;
handles.totRun = str2double(get(handles.totRun_editText,'String'));
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

handles.ROCarr = []; handles.PRarr = []; handles.tEarr = [];
handles.runLevl_statText.String = [num2str(0) '/' num2str(handles.totRun)]; pause(.001);
for c1 = 1:handles.totRun
    SDCOR(hObject,handles);
    handles = guidata(hObject);
    handles.ROCarr = [handles.ROCarr handles.ROC];
    handles.PRarr = [handles.PRarr handles.PR];
    handles.tEarr = [handles.tEarr handles.tElapsed];
    
    set(handles.tempROCPR_statText,'String',[num2str(handles.ROC,'%0.3f') ' / ' num2str(handles.PR,'%0.3f')]);
    set(handles.tempTime_statText,'String',num2str(handles.tElapsed,'%0.3f'));
    handles.runLevl_statText.String = [num2str(c1) '/' num2str(handles.totRun)]; pause(.001);
end
handles.ROCavg = mean(handles.ROCarr); handles.ROCstd = std(handles.ROCarr);
handles.PRavg = mean(handles.PRarr); handles.PRstd = std(handles.PRarr);
handles.tEavg = mean(handles.tEarr);

set(handles.ROCPRavg_statText,'String',[num2str(handles.ROCavg,'%0.3f') ' / ' num2str(handles.PRavg,'%0.3f')]);
set(handles.ROCPRstd_statText,'String',[num2str(handles.ROCstd,'%0.3f') ' / ' num2str(handles.PRstd,'%0.3f')]);
set(handles.runTime_statText,'String',num2str(handles.tEavg,'%0.3f'));
msgbox('Process was conducted successfully!','Success');

handles.startCond = 0; hOact(hObject,eventdata,handles);
PCMact(hObject,eventdata,handles);
    
guidata(hObject,handles);

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

function totRun_editText_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clearAxes_pushBtn_Callback(hObject,eventdata,handles)

cla(handles.axes1); legend(handles.axes1,'off');

handles.dsName_statText.String = '';
handles.progLevl_statText.String = '';
handles.tempROCPR_statText.String = '';
handles.tempTime_statText.String = '';
handles.runLevl_statText.String = '';
handles.ROCPRavg_statText.String = '';
handles.ROCPRstd_statText.String = '';
handles.runTime_statText.String = '';
handles.origK_val_statText.String = '';

guidata(hObject,handles);

function Hclear(hObject,eventdata,handles)

Harr = {'dsName','chunkSz','PCvarRat','alphaMemb','betaPrun','sampRate','totRun','PCM','PSO_particleNo','PSO_maxIter',...
    'PSO_W','PSO_C1','PSO_C2','PSO_alpha','manuEps','manuMnPt','epsCoeff','paramSampDS','paramCostArrSamp','origEps','origMnPt','origK',...
    'sampInd','idxSamp','mahalScores','idxFin','ROCarr','ROCavg','ROCstd','PRarr','PRavg','PRstd','tEarr','tEavg'};

for c1 = 1:numel(Harr)
    if isfield(handles,Harr{c1})
        handles = setfield(handles,Harr{c1},[]);
    end
end

guidata(hObject,handles);

function sampRate_editText_Callback(hObject,eventdata,handles)

function totRun_editText_Callback(hObject,eventdata,handles)

function PCvarRat_editText_Callback(hObject,eventdata,handles)

function PCvarRat_editText_KeyPressFcn(hObject,eventdata,handles)

function sampRate_editText_KeyPressFcn(hObject,eventdata,handles)

function ROCPRavg_statText_CreateFcn(hObject,eventdata,handles)

function chunkSz_editText_Callback(hObject,eventdata,handles)

function chunkSz_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function particleNo_editText_Callback(hObject,eventdata,handles)

function particleNo_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxIter_editText_Callback(hObject, eventdata, handles)

function maxIter_editText_CreateFcn(hObject, eventdata, handles)

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

function epsCoef_editText_Callback(hObject,eventdata,handles)

function epsCoef_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

function saveWork_pushBtn_Callback(hObject,eventdata,handles)

if isfield(handles,'ROCarr') && ~isempty(handles.ROCarr)
    Hsave = saveWork(hObject,eventdata,handles);
    uisave({'Hsave'},['..\results\','SDCOR(noVisDsk)_result_$',handles.dsName,'$_ROC=',num2str(handles.ROCavg,'%0.3f'),...
        '_ROCstd=',num2str(handles.ROCstd,'%0.3f'),'_PR=',num2str(handles.PRavg,'%0.3f'),'_PRstd=',num2str(handles.PRstd,'%0.3f'),...
        '_totRun=',num2str(handles.totRun),'_Time=',num2str(handles.tEavg,'%0.3f'),'.mat']);
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
    
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox('\fontsize{10}File was loaded successfully! Please load the dataset separately for a fresh test.','Success',CreateStruct));
end

function [Hsave] = saveWork(hObject,eventdata,handles)

Hsave = struct('dsName',handles.dsName,'chunkSz',handles.chunkSz,'PCvarRat',handles.PCvarRat*100,'alphaMemb',handles.alphaMemb,'betaPrun',handles.betaPrun,...
    'sampRate',handles.sampRate*100,'totRun',handles.totRun,'PCM',handles.PCM,'PSO_particleNo',handles.PSO_particleNo,...
    'PSO_maxIter',handles.PSO_maxIter,'PSO_W',handles.PSO_W,'PSO_C1',handles.PSO_C1,'PSO_C2',handles.PSO_C2,'PSO_alpha',handles.PSO_alpha,'manuEps',handles.manuEps,...
    'manuMnPt',handles.manuMnPt,'epsCoeff',handles.epsCoeff,'paramSampDS',handles.paramSampDS,'paramCostArrSamp',{handles.paramCostArrSamp},...
    'origEps',handles.origEps,'origMnPt',handles.origMnPt,'origK',handles.origK,'sampInd',handles.sampInd,'idxSamp',handles.idxSamp,'mahalScores',handles.mahalScores,...
    'idxFin',handles.idxFin,'ROCarr',handles.ROCarr,'ROCavg',handles.ROCavg,'ROCstd',handles.ROCstd,'PRarr',handles.PRarr,'PRavg',handles.PRavg,'PRstd',handles.PRstd,...
    'tEarr',handles.tEarr,'tEavg',handles.tEavg);

function [] = loadWork(hObject, eventdata, handles, Hsave)

clearAxes_pushBtn_Callback(hObject,eventdata,handles); % clearing workspace before loading the saved result
Hclear(hObject,eventdata,handles); handles = guidata(hObject);

handles.startCond = 0; hOact(hObject, eventdata, handles);
handles.dsName_statText.String = Hsave.dsName; handles.dsName = Hsave.dsName;

handles.labDS = [];
handles.chunkSz_editText.String = num2str(Hsave.chunkSz); handles.chunkSz = Hsave.chunkSz;
handles.PCvarRat_editText.String = num2str(Hsave.PCvarRat); handles.PCvarRat = Hsave.PCvarRat/100;
handles.alphaMemb_editText.String = num2str(Hsave.alphaMemb); handles.alphaMemb = Hsave.alphaMemb;
handles.betaPrun_editText.String = num2str(Hsave.betaPrun); handles.betaPrun = Hsave.betaPrun;
handles.sampRate_editText.String = num2str(Hsave.sampRate); handles.sampRate = Hsave.sampRate/100;
handles.totRun_editText.String = num2str(Hsave.totRun); handles.totRun = Hsave.totRun;

handles.PCM = Hsave.PCM;
switch Hsave.PCM
    case 'PSO_pcm_radioBtn'
        handles.PSO_pcm_radioBtn.Value = 1;
        
        axes(handles.axes1);
        plot(Hsave.paramCostArrSamp{1},'-r'); grid on;
        legend(sprintf('PSO costArr for SampDS\nEps=%0.3f, MinPts=%d',Hsave.paramSampDS(1),Hsave.paramSampDS(2)),'location','best');
        pause(.001);
        
    case 'manu_pcm_radioBtn'
        handles.manu_pcm_radioBtn.Value = 1;
        cla(handles.axes1);
        
end
PCMact(hObject,eventdata,handles);

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

handles.paramCostArrSamp = Hsave.paramCostArrSamp;
handles.paramSampDS = Hsave.paramSampDS;
handles.origEps = Hsave.origEps;
handles.origMnPt = Hsave.origMnPt;
handles.origK = Hsave.origK;
handles.sampInd = Hsave.sampInd;
handles.idxSamp = Hsave.idxSamp;

handles.mahalScores = Hsave.mahalScores;
handles.idxFin = Hsave.idxFin;
handles.ROCarr = Hsave.ROCarr; handles.ROCavg = Hsave.ROCavg; handles.ROCstd = Hsave.ROCstd;
handles.PRarr = Hsave.PRarr; handles.PRavg = Hsave.PRavg; handles.PRstd = Hsave.PRstd;
handles.tEarr = Hsave.tEarr; handles.tEavg = Hsave.tEavg;
handles.ROCPRavg_statText.String = [num2str(Hsave.ROCavg,'%0.3f') ' / ' num2str(Hsave.PRavg,'%0.3f')];
handles.ROCPRstd_statText.String = [num2str(Hsave.ROCstd,'%0.3f') ' / ' num2str(Hsave.PRstd,'%0.3f')];
handles.runTime_statText.String = num2str(Hsave.tEavg,'%0.3f');

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
    mainButns_Act(hObject,eventdata,handles,0);
    DBSCANparamCM_Act(hObject,eventdata,handles,0);
    
else
    SDCOR_InitParam_Act(hObject,eventdata,handles,1);
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
    handles.totRun_editText.Enable = 'off';
    handles.blckSzlim_editText.Enable = 'off';
    
else
    handles.chunkSz_editText.Enable = 'on';
    handles.PCvarRat_editText.Enable = 'on';
    handles.alphaMemb_editText.Enable = 'on';
    handles.betaPrun_editText.Enable = 'on';
    handles.sampRate_editText.Enable = 'on';
    handles.totRun_editText.Enable = 'on';
    handles.blckSzlim_editText.Enable = 'on';
    
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

function betaPrun_editText_Callback(hObject,eventdata,handles)

function betaPrun_editText_CreateFcn(hObject,eventdata,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nonUnifSamp_chckbx.
function nonUnifSamp_chckbx_Callback(hObject, eventdata, handles)
% hObject    handle to nonUnifSamp_chckbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonUnifSamp_chckbx
