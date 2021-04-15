%% Calculating accuracy and runtime values for DOLPHIN in multiple output case

% clear; close all; clc; warning off

% loading the labels file
[FileName,PathName] = uigetfile('*.labels','Please select a labels file','..\data and output\');
fileDIR = [PathName FileName]; [~,fName,~] = fileparts(FileName);
labels = load(fileDIR);

% loading the scores (outlier indices) files
[FileName,PathName] = uigetfile('Outliers_*.txt','Please select Outliers_* files','..\data and output\','MultiSelect','on');
if iscell(FileName); filCnt = numel(FileName); else; FileNameTemp{1} = FileName; FileName = FileNameTemp; filCnt = 1; end
OLarr = cell(filCnt,1);
for d1 = 1:filCnt
    fileDIR = [PathName FileName{d1}];
    OLarr{d1} = load(fileDIR);
end

OLs = cell2mat(OLarr);
OLsUnq = unique(cell2mat(OLarr));

scores = zeros(numel(labels),1);
for d2 = 1:numel(OLsUnq)
    scores(OLsUnq(d2)) = sum(OLs==OLsUnq(d2));
end

% loading the output files for calculation of the runtime
[FileName,PathName] = uigetfile('output_*.txt','Please select output_* files','..\data and output\','MultiSelect','on');
if iscell(FileName); filCnt = numel(FileName); else; FileNameTemp{1} = FileName; FileName = FileNameTemp; filCnt = 1; end
Tarr = zeros(filCnt,1);
for d1 = 1:filCnt
    fileDIR = [PathName FileName{d1}];
    fid = fopen(fileDIR,'r');
    while ~feof(fid)
        line = fgetl(fid);
        h1 = strfind(line,'Tempo elaborazione totale: ');
        if h1
            h2 = strfind(line,'sec');
            Tarr(d1) = str2double(line(28:h2-2));
            break
        end
    end
end
tElps = max(Tarr);

% computing ROC and PR accuracy values
[~,~,~,ROC] = perfcurve(labels,scores,1);
[~,~,~,PR] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('DOLPHIN results for %s:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',fName,ROC,PR,tElps);
save(['res_DOLPHIN_' fName '.mat'],'ROC','PR','tElps');

