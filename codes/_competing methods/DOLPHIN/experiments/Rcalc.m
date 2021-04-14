clear; close all; clc; warning off


%% processing real-life datasets

Eps = .01; delta = .1; alpha = [1:10]./100; sigma = .01;
alphaK = numel(alpha);

addpath('G:\Dropbox\Researches\Implementations\SDCOR\datasets\realDS')

% Mammography
load('mammography_(11183by6_260o).mat');
R_mammography = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_mammography(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Adult
load('adult_(38323by6_1168o).mat');
R_adult = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_adult(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Shuttle
load('shuttle_(49097by9_3511o).mat');
R_shuttle = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_shuttle(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Smtp
load('smtp_(95156by3_30o).mat');
R_smtp = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_smtp(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Skin
load('skin_(199283by3_5085o).mat');
R_skin = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_skin(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% CreditCardFraud
load('creditCardFraud_(284807by29_492o).mat');
R_creditCardFraud = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_creditCardFraud(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% ForestCover
load('forestCover_(286048by10_2747o).mat');
R_forestCover = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_forestCover(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Http
load('http_(567498by3_2211o).mat');
R_http = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_http(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% Hepc
load('hepc_(2003171by7_5123o).mat');
R_hepc = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_hepc(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% save(['R_realData_(n=' num2str(n) ',q=' num2str(q) ',a=' num2str(a) ')']);

%% processing synthetic datasets

%** considering "Accuracy and Stability Analysis" datasets

Eps = .01; delta = .1; alpha = [1:10]./100; sigma = .01;
alphaK = numel(alpha);

addpath('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\datasets\synthDS\Accuracy and Stability Analysis')

% accrData1
load('accrData1_(500000by30_5000o).mat');
R_accrData1 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_accrData1(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% accrData2
load('accrData2_(1000000by40_10000o).mat');
R_accrData2 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_accrData2(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% accrData3
load('accrData3_(1500000by50_15000o).mat');
R_accrData3 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_accrData3(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% accrData4
load('accrData4_(2000000by60_20000o).mat');
R_accrData4 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_accrData4(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% save(['R_accrData_(n=' num2str(n) ',q=' num2str(q) ',a=' num2str(a) ')']);

%** considering "Tolerance to a High Number of Outliers" datasets

Eps = .01; delta = .1; alpha = [5:5:50]./100; sigma = .01;
alphaK = numel(alpha);

addpath('G:\Dropbox Aux. Folder\Researches\Implementations\SDCOR\datasets\synthDS\Tolerance to a High Number of Outliers\')

% tolrData01
load('tolrData01_(30000by2_10000o).mat');
R_tolrData01 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData01(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData02
load('tolrData02_(32000by2_12000o).mat');
R_tolrData02 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData02(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData03
load('tolrData03_(34000by2_14000o).mat');
R_tolrData03 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData03(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData04
load('tolrData04_(36000by2_16000o).mat');
R_tolrData04 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData04(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData05
load('tolrData05_(38000by2_18000o).mat');
R_tolrData05 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData05(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData06
load('tolrData06_(40000by2_20000o).mat');
R_tolrData06 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData06(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData07
load('tolrData07_(42000by2_22000o).mat');
R_tolrData07 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData07(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData08
load('tolrData08_(44000by2_24000o).mat');
R_tolrData08 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData08(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData09
load('tolrData09_(46000by2_26000o).mat');
R_tolrData09 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData09(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData10
load('tolrData10_(48000by2_28000o).mat');
R_tolrData10 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData10(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% tolrData11
load('tolrData11_(50000by2_30000o).mat');
R_tolrData11 = zeros(1,alphaK);
for c1 = 1:alphaK
    [R_tolrData11(c1)] = DolphinParamEstim(X,Eps,delta,alpha(c1),sigma);
end

% save(['R_tolrData_(n=' num2str(n) ',q=' num2str(q) ',a=' num2str(a) ')']);


