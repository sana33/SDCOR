%% Calculating accuracy and runtime values for DOLPHIN in single output case

% clear; close all; clc; warning off

% loading the labels file
[FileName,PathName] = uigetfile('*.labels','Please select a labels file','..\data and output\');
fileDIR = [PathName FileName]; [~,fName,~] = fileparts(FileName);
labels = load(fileDIR);

% loading the scores (outliers) file
[FileName,PathName] = uigetfile('Outliers_*.txt','Please select Outliers_* files','..\data and output\');
fileDIR = [PathName FileName]; scorTemp = load(fileDIR); scores = zeros(numel(labels),1); scores(scorTemp) = 1;

% computing ROC and PR accuracy values
[~,~,~,ROC] = perfcurve(labels,scores,1);
[~,~,~,PR] = perfcurve(labels,scores,1,'XCrit','reca','YCrit','prec');

fprintf('DOLPHIN results for %s:\t\tROC = %0.3f\t\tPR = %0.3f\n\n',fName,ROC,PR);

