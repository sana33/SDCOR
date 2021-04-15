%% Calculating accuracy values for DOLPHIN in single output case

% clear; close all; clc; warning off

% loading the labels file
[FileName,PathName] = uigetfile('*.labels','Please select a labels file','..\data and output\');
fileDIR = [PathName FileName]; [~,fName,~] = fileparts(FileName);
labels = load(fileDIR);

% loading the scores (outlier indices) file
[FileName,PathName] = uigetfile('Outliers_*.txt','Please select Outliers_* files','..\data and output\');
fileDIR = [PathName FileName]; scorTemp = load(fileDIR); scores = zeros(numel(labels),1); scores(scorTemp) = 1;

% loading the output file for calculation of the runtime
[FileName,PathName] = uigetfile('output_*.txt','Please select output_* files','..\data and output\');
fileDIR = [PathName FileName];
fid = fopen(fileDIR,'r');
while ~feof(fid)
	line = fgetl(fid);
	h1 = strfind(line,'Tempo elaborazione totale: ');
	if h1
		h2 = strfind(line,'sec');
		tElps = str2double(line(28:h2-2));
		break
	end
end

% computing ROC and PR accuracy values
[~,~,~,ROC] = perfcurve(labels,scores,1);
[~,~,~,PR] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('DOLPHIN results for %s:\t\tROC = %0.3f\t\tPR = %0.3f\t\telpsTime = %0.3f sec\n\n',fName,ROC,PR,tElps);
save(['res_DOLPHIN_' fName '.mat'],'ROC','PR','tElps');
