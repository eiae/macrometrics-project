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
addpath('functions_dfm')
addpath('functions_var')

% seeds
seed = 1234;  
rng(seed);


%% data 

% load data
datatable = readtable('data.xlsx', ReadVariableNames=true);  
names = datatable.Properties.VariableNames;
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

select = [1 2 3 4 5 7];  % select specific variables
T = length(dates);
N = length(select);


%% preprocess data
run('preprocess')


%% run models
run('dfm');
run('var');


%% plot charts

