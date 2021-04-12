
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

function MAIN_OpeningFcn(hO, eventdata, H, varargin)

H.output = hO;

warning off;
PCMact(hO,eventdata,H);

% Update H structure
guidata(hO, H);

function varargout = MAIN_OutputFcn(hO,eventdata,H) 

varargout{1} = H.output;

function load_pushBtn_Callback(hO,eventdata,H)

clearAxes_pushBtn_Callback(hO,eventdata,H); % clearing workspace before loading new data

[FileName,PathName] = uigetfile('*.mat', 'Select the dataset along with outlier labels, all as a single MAT-file','..\datasets\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    H.labDS = matfile([PathName FileName]);
    [H.n,H.p] = size(H.labDS,'X');
    
    [~,H.dsName,~] = fileparts(FileName);
    H.dsName_statText.String = H.dsName;
    H.chunkSz_editText.String = ceil(.1*H.n);
    
    H.manu_pcm_radioBtn.Value = 1; PCMact(hO,eventdata,H);
    msgbox('File was loaded successfully!','Success');
end

guidata(hO,H);

function start_pushBtn_Callback(hO,eventdata,H)

global BLK_SZ_LIM

%------- error handing -------%
if ~isfield(H,'labDS') || isempty(H.labDS)
    errordlg('Dataset file not found! Please load the input data first!','File Error');
    return
end
%-----------------------------%

H.startCond = 1; hOact(hO,eventdata,H);

H.dsName_statText.String = H.dsName;
H.chunkSz = str2double(get(H.chunkSz_editText,'String'));
H.PCvarRat = str2double(get(H.PCvarRat_editText,'String'))/100;
H.alphaMemb = str2double(get(H.alphaMemb_editText,'String'));
H.betaPrun = str2double(get(H.betaPrun_editText,'String'));
H.sampRate = str2double(get(H.sampRate_editText,'String'))/100;
H.maxRun = str2double(get(H.maxRun_editText,'String'));
BLK_SZ_LIM = str2double(get(H.blckSzlim_editText,'String'));

H.PCM = get(get(H.PCM_radioBtnGroup,'SelectedObject'),'tag');
H.PSO_particleNo = str2double(get(H.particleNo_editText,'String'));
H.PSO_maxIter = str2double(get(H.maxIter_editText,'String'));
H.PSO_W = str2double(get(H.W_editText,'String'));
H.PSO_C1 = str2double(get(H.C1_editText,'String'));
H.PSO_C2 = str2double(get(H.C2_editText,'String'));
H.PSO_alpha = str2double(get(H.alpha_editText,'String'));
H.manuEps = str2double(get(H.manuEps_editText,'String'));
if isnan(H.manuEps); H.manuEps = 0; end
H.manuMnPt = str2double(get(H.manuMnPt_editText,'String'));
if isnan(H.manuMnPt); H.manuMnPt = floor(log(H.n)); H.manuMnPt_editText.String = num2str(H.manuMnPt); end
H.epsCoeff = str2double(get(H.epsCoef_editText,'String'));

H.ROCarr = []; H.PRarr = []; H.tEarr = [];
H.runLevl_statText.String = [num2str(0) '/' num2str(H.maxRun)]; pause(.001);
for c1 = 1:H.maxRun
    SDCOR(hO,H);
    H = guidata(hO);
    H.ROCarr = [H.ROCarr H.ROC];
    H.PRarr = [H.PRarr H.PR];
    H.tEarr = [H.tEarr H.tElapsed];
    
    set(H.tempROCPR_statText,'String',[num2str(H.ROC,'%0.3f') ' / ' num2str(H.PR,'%0.3f')]);
    set(H.tempTime_statText,'String',num2str(H.tElapsed,'%0.3f'));
    H.runLevl_statText.String = [num2str(c1) '/' num2str(H.maxRun)]; pause(.001);
end
H.ROCavg = mean(H.ROCarr); H.ROCstd = std(H.ROCarr);
H.PRavg = mean(H.PRarr); H.PRstd = std(H.PRarr);
H.tEavg = mean(H.tEarr);

set(H.ROCPRavg_statText,'String',[num2str(H.ROCavg,'%0.3f') ' / ' num2str(H.PRavg,'%0.3f')]);
set(H.ROCPRstd_statText,'String',[num2str(H.ROCstd,'%0.3f') ' / ' num2str(H.PRstd,'%0.3f')]);
set(H.runTime_statText,'String',num2str(H.tEavg,'%0.3f'));
msgbox('Process was conducted successfully!','Success');

H.startCond = 0; hOact(hO,eventdata,H);
PCMact(hO,eventdata,H);
    
guidata(hO,H);

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

function maxRun_editText_CreateFcn(hO, eventdata, H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function clearAxes_pushBtn_Callback(hO,eventdata,H)

cla(H.axes1); legend(H.axes1,'off');

H.dsName_statText.String = '';
H.progLevl_statText.String = '';
H.tempROCPR_statText.String = '';
H.tempTime_statText.String = '';
H.runLevl_statText.String = '';
H.ROCPRavg_statText.String = '';
H.ROCPRstd_statText.String = '';
H.runTime_statText.String = '';
H.origK_val_statText.String = '';

guidata(hO,H);

function sampRate_editText_Callback(hO,eventdata,H)

function maxRun_editText_Callback(hO,eventdata,H)

function PCvarRat_editText_Callback(hO,eventdata,H)

function PCvarRat_editText_KeyPressFcn(hO,eventdata,H)

function sampRate_editText_KeyPressFcn(hO,eventdata,H)

function ROCPRavg_statText_CreateFcn(hO,eventdata,H)

function chunkSz_editText_Callback(hO,eventdata,H)

function chunkSz_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function particleNo_editText_Callback(hO,eventdata,H)

function particleNo_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function maxIter_editText_Callback(hO, eventdata, H)

function maxIter_editText_CreateFcn(hO, eventdata, H)

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

function epsCoef_editText_Callback(hO,eventdata,H)

function epsCoef_editText_CreateFcn(hO,eventdata,H)

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

function saveWork_pushBtn_Callback(hO,eventdata,H)

if isfield(H,'ROCarr')
    Hsave = saveWork(hO,eventdata,H);
    uisave({'Hsave'},['..\results\','SDCOR(noVisDsk)_result_$',H.dsName,'$_ROC=',num2str(H.ROCavg,'%0.3f'),...
        '_ROCstd=',num2str(H.ROCstd,'%0.3f'),'_PR=',num2str(H.PRavg,'%0.3f'),'_PRstd=',num2str(H.PRstd,'%0.3f'),...
        '_maxRun=',num2str(H.maxRun),'_Time=',num2str(H.tEavg,'%0.3f'),'.mat']);
    msgbox('File was saved successfully!','Success');
else
    msgbox('Sorry! There is not any clear run info to be saved!','Failure','error');
end

function loadWork_pushBtn_Callback(hO,eventdata,H)

clearAxes_pushBtn_Callback(hO,eventdata,H); % clearing workspace before loading the saved result

[FileName,PathName] = uigetfile('*.mat', 'Select the saved workspace to be loaded','..\results\');
if ~FileName
    msgbox('Sorry! No file was loaded!','Failure','error');
else
    Hsave = importdata([PathName FileName]);
    loadWork(hO, eventdata, H, Hsave);
    
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox('\fontsize{10}File was loaded successfully! Please load the dataset separately for a fresh test.','Success',CreateStruct));
end

function [Hsave] = saveWork(hO,eventdata,H)

global BLK_SZ_LIM

Hsave = struct('dsName',H.dsName,'chunkSz',H.chunkSz,'PCvarRat',H.PCvarRat*100,'alphaMemb',H.alphaMemb,'betaPrun',H.betaPrun,...
    'sampRate',H.sampRate*100,'maxRun',H.maxRun,'BLK_SZ_LIM',BLK_SZ_LIM,'PCM',H.PCM,'PSO_particleNo',H.PSO_particleNo,...
    'PSO_maxIter',H.PSO_maxIter,'PSO_W',H.PSO_W,'PSO_C1',H.PSO_C1,'PSO_C2',H.PSO_C2,'PSO_alpha',H.PSO_alpha,'manuEps',H.manuEps,...
    'manuMnPt',H.manuMnPt,'epsCoeff',H.epsCoeff,'paramSampDS',H.paramSampDS,'paramCostArrSamp',{H.paramCostArrSamp},...
    'origEps',H.origEps,'origMnPt',H.origMnPt,'origK',H.origK,'sampInd',H.sampInd,'idxSamp',H.idxSamp,'mahalScores',H.mahalScores,...
    'idxFin',H.idxFin,'ROCarr',H.ROCarr,'ROCavg',H.ROCavg,'ROCstd',H.ROCstd,'PRarr',H.PRarr,'PRavg',H.PRavg,'PRstd',H.PRstd,...
    'tEarr',H.tEarr,'tEavg',H.tEavg);

function [] = loadWork(hO, eventdata, H, Hsave)

global BLK_SZ_LIM

H.startCond = 0; hOact(hO, eventdata, H);
H.dsName_statText.String = Hsave.dsName;

H.labDS = [];
H.chunkSz_editText.String = num2str(Hsave.chunkSz);
H.PCvarRat_editText.String = num2str(Hsave.PCvarRat);
H.alphaMemb_editText.String = num2str(Hsave.alphaMemb);
H.betaPrun_editText.String = num2str(Hsave.betaPrun);
H.sampRate_editText.String = num2str(Hsave.sampRate);
H.maxRun_editText.String = num2str(Hsave.maxRun);
BLK_SZ_LIM = Hsave.BLK_SZ_LIM; H.blckSzlim_editText.String = num2str(BLK_SZ_LIM);

switch Hsave.PCM
    case 'PSO_pcm_radioBtn'
        H.PSO_pcm_radioBtn.Value = 1;
        
        axes(H.axes1);
        plot(Hsave.paramCostArrSamp{1},'-r'); grid on;
        legend(sprintf('PSO costArr for SampDS\nEps=%0.3f, MinPts=%d',Hsave.paramSampDS(1),Hsave.paramSampDS(2)),'location','best');
        pause(.001);
        
    case 'manu_pcm_radioBtn'
        H.manu_pcm_radioBtn.Value = 1;
        cla(H.axes1);
        
end
PCMact(hO,eventdata,H);

H.particleNo_editText.String = num2str(Hsave.PSO_particleNo);
H.maxIter_editText.String = num2str(Hsave.PSO_maxIter);
H.W_editText.String = num2str(Hsave.PSO_W);
H.C1_editText.String = num2str(Hsave.PSO_C1);
H.C2_editText.String = num2str(Hsave.PSO_C2);
H.alpha_editText.String = num2str(Hsave.PSO_alpha);
H.manuEps_editText.String = num2str(Hsave.manuEps);
H.manuMnPt_editText.String = num2str(Hsave.manuMnPt);
H.epsCoef_editText.String = num2str(Hsave.epsCoeff);

H.paramCostArrSamp = Hsave.paramCostArrSamp;
H.paramSampDS = Hsave.paramSampDS;
H.origEps = Hsave.origEps;
H.origMnPt = Hsave.origMnPt;
H.origK = Hsave.origK;
H.sampInd = Hsave.sampInd;
H.idxSamp = Hsave.idxSamp;

H.mahalScores = Hsave.mahalScores;
H.idxFin = Hsave.idxFin;
H.ROCarr = Hsave.ROCarr; H.ROCavg = Hsave.ROCavg; H.ROCstd = Hsave.ROCstd;
H.PRarr = Hsave.PRarr; H.PRavg = Hsave.PRavg; H.PRstd = Hsave.PRstd;
H.tEarr = Hsave.tEarr; H.tEavg = Hsave.tEavg;
H.ROCPRavg_statText.String = [num2str(H.ROCavg,'%0.3f') ' / ' num2str(H.PRavg,'%0.3f')];
H.ROCPRstd_statText.String = [num2str(H.ROCstd,'%0.3f') ' / ' num2str(H.PRstd,'%0.3f')];
H.runTime_statText.String = num2str(H.tEavg,'%0.3f');
H.origK_val_statText.String = num2str(Hsave.origK);

guidata(hO,H);

function PCMact(hO,eventdata,H)

switch H.PCM_radioBtnGroup.SelectedObject.Tag
    case 'Kgraph_pcm_radioBtn'
        H.makeManu_pushBtn.Enable =  'off';
        
        H.particleNo_editText.Enable = 'off';
        H.maxIter_editText.Enable = 'off';
        H.W_editText.Enable = 'off';
        H.C1_editText.Enable = 'off';
        H.C2_editText.Enable = 'off';
        H.alpha_editText.Enable = 'off';
        H.manuEps_editText.Enable = 'off';
        H.manuMnPt_editText.Enable = 'on';
        H.epsCoef_editText.Enable = 'off';
        
    case 'PSO_pcm_radioBtn'
        H.makeManu_pushBtn.Enable =  'on';
        
        H.particleNo_editText.Enable = 'on';
        H.maxIter_editText.Enable = 'on';
        H.W_editText.Enable = 'on';
        H.C1_editText.Enable = 'on';
        H.C2_editText.Enable = 'on';
        H.alpha_editText.Enable = 'on';
        H.manuEps_editText.Enable = 'off';
        H.manuMnPt_editText.Enable = 'on';
        H.epsCoef_editText.Enable = 'on';
        
    case 'manu_pcm_radioBtn'
        H.makeManu_pushBtn.Enable =  'off';
        
        H.particleNo_editText.Enable = 'off';
        H.maxIter_editText.Enable = 'off';
        H.W_editText.Enable = 'off';
        H.C1_editText.Enable = 'off';
        H.C2_editText.Enable = 'off';
        H.alpha_editText.Enable = 'off';
        H.manuEps_editText.Enable = 'on';
        H.manuMnPt_editText.Enable = 'on';
        H.epsCoef_editText.Enable = 'on';
        
end

function hOact(hO, eventdata, H)

if H.startCond
    SDCOR_InitParam_Act(hO,eventdata,H,0);
    clearAxes_pushBtn_Callback(hO,eventdata,H);
    mainButns_Act(hO,eventdata,H,0);
    DBSCANparamCM_Act(hO,eventdata,H,0);
    
else
    SDCOR_InitParam_Act(hO,eventdata,H,1);
    DBSCANparamCM_Act(hO,eventdata,H,1);
    mainButns_Act(hO,eventdata,H,1);
    
end

function SDCOR_InitParam_Act(hO, eventdata, H, actCond)

if ~actCond
    H.chunkSz_editText.Enable = 'off';
    H.PCvarRat_editText.Enable = 'off';
    H.alphaMemb_editText.Enable = 'off';
    H.betaPrun_editText.Enable = 'off';
    H.sampRate_editText.Enable = 'off';
    H.maxRun_editText.Enable = 'off';
    H.blckSzlim_editText.Enable = 'off';
    
else
    H.chunkSz_editText.Enable = 'on';
    H.PCvarRat_editText.Enable = 'on';
    H.alphaMemb_editText.Enable = 'on';
    H.betaPrun_editText.Enable = 'on';
    H.sampRate_editText.Enable = 'on';
    H.maxRun_editText.Enable = 'on';
    H.blckSzlim_editText.Enable = 'on';
    
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
    H.Kgraph_pcm_radioBtn.Enable = 'off';
    H.PSO_pcm_radioBtn.Enable = 'off';
    H.manu_pcm_radioBtn.Enable = 'off';
    H.makeManu_pushBtn.Enable =  'off';
    
    H.particleNo_editText.Enable = 'off';
    H.maxIter_editText.Enable = 'off';
    H.W_editText.Enable = 'off';
    H.C1_editText.Enable = 'off';
    H.C2_editText.Enable = 'off';
    H.alpha_editText.Enable = 'off';
    H.manuEps_editText.Enable = 'off';
    H.manuMnPt_editText.Enable = 'off';
    H.epsCoef_editText.Enable = 'off';
else
    H.Kgraph_pcm_radioBtn.Enable = 'on';
    H.PSO_pcm_radioBtn.Enable = 'on';
    H.manu_pcm_radioBtn.Enable = 'on';
    H.makeManu_pushBtn.Enable =  'on';
    
    H.particleNo_editText.Enable = 'on';
    H.maxIter_editText.Enable = 'on';
    H.W_editText.Enable = 'on';
    H.C1_editText.Enable = 'on';
    H.C2_editText.Enable = 'on';
    H.alpha_editText.Enable = 'on';
    H.manuEps_editText.Enable = 'on';
    H.manuMnPt_editText.Enable = 'on';
    H.epsCoef_editText.Enable = 'on';
end

function Kgraph_pcm_radioBtn_Callback(hO, eventdata, H)

H.manuEps_editText.String = '';

if isfield(H,'p')
    H.manuMnPt_editText.String = 10*H.p;
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    msgCont = '\fontsize{10} To avoid outliers, we set the {\it{MinPts}} quantity to {\bf{10{\cdot}p}}. You can change it at your will!';
    uiwait(msgbox(msgCont,'Alert!','Help',CreateStruct));
end

PCMact(hO,eventdata,H);
   
function PSO_pcm_radioBtn_Callback(hO,eventdata,H)

H.manuEps_editText.String = '';

if isfield(H,'p')
    H.manuMnPt_editText.String = 10*H.p;
    CreateStruct.Interpreter = 'tex'; CreateStruct.WindowStyle = 'modal';
    msgCont = ['\fontsize{10} PSO searches in predetermined ranges by the user for finding the optimal {\it{Eps}} and {\it{MinPts}} ',...
        'quantities. The {\it{Eps}} range and the lower bound for the {\it{MinPts}} range are defined automatically; for the {\it{MinPts}} ',...
        'upper bound, we set it to {\bf{10{\cdot}p}} for your convenience.\newline\newline {\color{red}\bf{Note:}} For being more prudent ',...
        'in the case of high noisy datasets, you can change it to greater values at your will! But please do NOT go much further as for ',...
        'vary large {\it{MinPts}} values, the following {\it{k}}-dist graph tends to become horizontal and without any distinguishable ',...
        'valleys; therefore, it might not be possible to choose the right value for the {\it{MinPts}} parameter.'];
    uiwait(msgbox(msgCont,'Alert!','Help',CreateStruct));
end

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

H.startCond = 0; hOact(hO, eventdata, H);
PCMact(hO,eventdata,H);

guidata(hO,H);

function blckSzlim_editText_Callback(hO,eventdata,H)

function blckSzlim_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end

function betaPrun_editText_Callback(hO,eventdata,H)

function betaPrun_editText_CreateFcn(hO,eventdata,H)

if ispc && isequal(get(hO,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hO,'BackgroundColor','white');
end
