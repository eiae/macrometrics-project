%% Macroeconometrics project
% -------------------------------------------------------------------------
% Forecasting performance between models and relative to ECB B/MPE
% - Mixed-frequency Bayesian Dynamic Factor Model
% - Mixed-frequency Bayesian Vector Autoregressive Model
% @ Erik Andres Escayola 
% -------------------------------------------------------------------------


%% housekeepiing
clc 
clear 
close all

% paths
tmp = matlab.desktop.editor.getActive;  
cwd = fileparts(tmp.Filename);
cd(cwd);  
addpath('functions')
addpath('functions_dfm')
addpath('functions_var')

% seeds
seed = 1234;  
rng(seed);


%% data 

% load data
datatable = readtable('data.xlsx', ReadVariableNames=true);  
names = datatable.Properties.VariableNames(2:end);
dates = table2array(datatable(:,1));
data = table2array(datatable(:,2:end));

% variable description (Mnemonic -> Variable - Sector - Type - Frequency)
% 1  GDP -> Real GDP - Activity - Hard - Q	
% 2  PMI -> Composite output PMI - Activity - Soft - M
% 3  ESI ->	Economic sentiment - Activity - Soft - M
% 4  IP  ->	Industrial production manufacturing - Activity - Hard - M
% 5  RS	 -> Retail sales - Activity - Hard - M
% 6  EMP ->	Unemployment rate - Activity - Hard -M
% 7  EEXP->	Extra-euro area exports - Activity - Hard - M
% 8  CPI ->	Consumer prices - Prices - Hard - M
% 9  REER->	Real effective exchange rate - Prices - Hard - M
% 10 COM ->	Commodity spot aggregate - Prices - Hard - M
% 11 CISS->	Sovereign stress - Financial - Hard/Soft - M
% 12 EQTY->	EURO STOXX 50 - Financial - Hard - M
% 13 INT ->	Short-term interest rate & shadow rate - Financial - Hard/Soft - M
% 14 EPU ->	Economic policy uncertainty - Other - Soft - M
% 15 SBI -> Supply-bottleneck - Other - Soft - M
% Note: always put quarterly variables first

selected = [1 2 3 4 5 7];  % select specific variables
%selected = [1 2 3 4 5 7 8 9 10 11 12 14 15];
%selected = [1 4 5 6 7 9 10];
%selected = 1:size(data,2);
names = names(selected);

T = length(dates);
N = length(selected);

Q = 1;  % number of quarterly variables
M = N-Q;  % number of monthly variables

yraw = data(:,selected);  % observables
idx = 1-isnan(yraw);  % indicator to select filled values

prct = [10 90];  % percentiles for credible bands
H = 36;  % forecast horizon


%% preprocess data
yDFM = preprocessDFM(yraw,dates,names,selected,T,N,Q,M);forecast
yVAR = preprocessVAR(yraw,dates,names,selected,T,N,Q,M);

% no covid
% yraw = yraw(1:end-12*4,:);  
% T = size(yraw,1);
% idx = 1-isnan(yraw);
% dates = dates(1:end-12*4,:);
% y = preprocess(yraw,dates,names,selected,T,N,Q,M);


%% run models
run('DFM');
run('VAR');


%% ECB B/MPE projections

% load ECB data
datatableECB = readtable('ECB_projections.xlsx', ReadVariableNames=true);  
namesECB = datatableECB.Properties.VariableNames(2:end);
datesECB = table2array(datatableECB(:,1));
dataECB = table2array(datatableECB(:,2:end));

% plot latest data along projections
periodsBack = 5*4;  % years times quarters

figure;
plot(datesECB(end-periodsBack:end), dataECB(end-periodsBack:end,2), Color=[1 0.71 0], LineWidth=1.5)
hold on
plot(datesECB(end-periodsBack:end), dataECB(end-periodsBack:end,1), Color=[0 0.22 0.6], LineWidth=1.5)
axis tight
grid on 
legend(namesECB{2}, namesECB{1})

sgt = sgtitle('ECB Macroeconomic Projections on euro area real GDP growth', 'Interpreter','latex');
sgt.FontSize = 20;


%% comparison analysis

datesECBDFM = datesECB(2:end);  % match ECB and model dates due to lag for computing growth rates
dataECBDFM = dataECB(2:end,:);  % match ECB and model data due to lag for computing growth rates

% compare DFM and ECB
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
plot(datesECBDFM(end-periodsBack:end), dataECBDFM(end-periodsBack:end,2), Color=[1 0.71 0], LineWidth=1.5)
hold on
plot(datesECBDFM(end-periodsBack:end), dataECBDFM(end-periodsBack:end,1), Color=[0 0.22 0.6], LineWidth=1.5)
ylabel('ECB')

yyaxis right
plot(datesECBDFM(end-periodsBack:end), forecastQ(end-periodsBack:end), Color='k', LineWidth=1.5)
hFill = [datesECBDFM(end-periodsBack:end)' fliplr(datesECBDFM(end-periodsBack:end)')];
inBetween = [forecast_bandsQ(end-periodsBack:end,1)', fliplr(forecast_bandsQ(end-periodsBack:end,end)')];
fill(hFill , inBetween, 'r', FaceAlpha=0.2, LineStyle='none');
ylabel('DFM')

axis tight
grid on 
legend(namesECB{2}, namesECB{1}, 'Median forecast DFM', 'Credible bands 90%', Location='best')
sgt = sgtitle('Comparison ECB projections and DFM forecast', 'Interpreter','latex');
sgt.FontSize = 20;


% compare VAR to ECB



% compare DFM and VAR to ECB

% compute deltas and correlations
horizonECB = 4*3;  % 3-year forecast horizon of ECB 
forecastStart = length(datesECB) - length(datesECBDFM(end-horizonECB:end));

forecastDeltaECBMar24vsDFM = dataECBDFM(forecastStart:end,1) - forecastQ(forecastStart:end);
[fECBMar24vsDFM, xiECBMar24vsDFM] = ksdensity(forecastDeltaECBMar24vsDFM);

forecastCorrECBMar24vsDFM = corr(dataECBDFM(forecastStart:end,1), forecastQ(forecastStart:end));
fprintf('Forecast correlation between ECB Mar23 and DFM: %4.2f', forecastCorrECBMar24vsDFM)

figure
plot(xiECBMar24vsDFM,fECBMar24vsDFM, Color=[0 0.22 0.6], LineWidth=1.5);
title('Kernel density of forecast differences between ECB Mar23 MPE and DFM')

