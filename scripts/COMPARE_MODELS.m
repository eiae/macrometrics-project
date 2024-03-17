%% Comparison analysis between models

% Disclaimer: code written
% by Erik Andres Escayola

% clear
clc


%% forecast comparison exercise

% specs
forecastPeriodsBack = 12;  % months, minimum 12 months = 1 year
forecastPeriodsBackQ = forecastPeriodsBack/3;
H = 3;  % redefine horizon for exercie; months, minimum 3 months = 1 quarter
HQ = H/3;

% target variable with full sample to compare forecast and actual values
Tfull = size(yraw,1);  
datesfull = table2array(datatable(:,1));

% forecast comparison
minWindow = Tfull-forecastPeriodsBack; 
forecastPeriods = length(minWindow:3:Tfull-H);

RmseDFM = zeros(forecastPeriods, 1);
RmseVAR = zeros(forecastPeriods, 1);


%% forecast comparison between DFM and actual

% preprocess data and get actual with full sample
yDFMFull = preprocessDFM(yraw,selected,Tfull,N,Q);
obsFull = yDFMFull(:,1);
obsFullQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= Tfull  
    obsFullQ = [obsFullQ; obsFull(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% estimate model in expanding window 
counter = 1;
for window = minWindow:3:Tfull-H
    
    fprintf('Forecast iteration progress: %d%%\n',  round(counter/forecastPeriods*100,1))

    % select subsample
    subsample = yraw(1:window,:);
    T = size(subsample,1);
    idx = 1-isnan(subsample);
    dates = datesfull(1:window,:);

    % transformations model
    yDFM = preprocessDFM(subsample,selected,T,N,Q);

    % run model
    run(fullfile(workpath,'DFM.m'))
    %run(fullfile(workpath,'DFM_CHART.m'))  % uncomment to check estimation results in each forecast

    % compare forecast with actual data
    %actualValueDFM = obsFullQ(end-forecastPeriodsBackQ+1:end-forecastPeriodsBackQ+HQ);
    backShift = length(forecastQ(1:end-HQ));
    actualValueDFM = obsFullQ(backShift+1:backShift+HQ);
    forecastValueDFM = forecastQ(end-HQ+1:end);  %end-HQ+1:end 
    
    disp('actual vs forecast values for horizon')
    disp([actualValueDFM,forecastValueDFM])

    % compute rmse
    %error = actual_value - forecast_value;
    %RMSE= sqrt(mean(error.^2))
    RmseValueDFM = rmse(actualValueDFM, forecastValueDFM);

    
    % store rmse value
    RmseDFM(counter) = RmseValueDFM;

    counter = counter + 1;
end

% compute kernel density
[fRmseDFM, xiRmseDFM] = ksdensity(RmseDFM);


%% forecast comparison between DFM and actual

% preprocess data and get actual with full sample
yVARFull = preprocessVAR(yraw,selected,Tfull,N,Q);
obsAltFull = [NaN(L,1); 100*( yVARFull(L+1:end,1) - yVARFull(1:end-L,1) )];  % quarter-on-quarter growth rates based on monthly series of quartely data
obsAltFullQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= Tfull  
    obsAltFullQ = [obsAltFullQ; obsAltFull(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% estimate model in expanding window 
counter = 1;
for window = minWindow:3:Tfull-H
    
    fprintf('Forecast iteration progress: %d%%\n',  round(counter/forecastPeriods*100,1))

    % select subsample
    subsample = yraw(1:window,:);
    T = size(subsample,1);
    idx = 1-isnan(subsample);
    dates = datesfull(1:window,:);

    % transformations model
    yVAR = preprocessVAR(subsample,selected,T,N,Q);

    % run model
    run(fullfile(workpath,'VAR.m'))
    %run(fullfile(workpath,'VAR_CHART.m'))  % uncomment to check estimation results in each forecast

    % compare forecast with actual data
    %actualValueVAR = obsAltFullQ(end-forecastPeriodsBackQ+1:end-forecastPeriodsBackQ+HQ);
    backShift = length(forecastAltQ(1:end-HQ));
    actualValueVAR = obsAltFullQ(backShift+1:backShift+HQ);
    forecastValueVAR = forecastAltQ(end-HQ+1:end);

    disp('actual vs forecast values for horizon')
    disp([actualValueVAR,forecastValueVAR])
    
    % compute rmse
    %error = actual_value - forecast_value;
    %RMSE= sqrt(mean(error.^2))
    RmseValueVAR = rmse(actualValueVAR, forecastValueVAR);
    
    % store rmse value
    RmseVAR(counter) = RmseValueVAR;

    counter = counter + 1;
end

% compute kernel density
[fRmseVAR, xiRmseVAR] = ksdensity(RmseVAR);


%% charts

% RMSE values along forecast window
figure;
subplot(2,1,1)
plot(RmseDFM, Color=colorDFM, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
title('RMSE DFM')

subplot(2,1,2)
plot(RmseVAR, Color=colorVAR, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
title('RMSE VAR')

% RMSE kernel densities
figure;
subplot(2,1,1)
plot(xiRmseDFM, fRmseDFM, Color=colorDFM, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
title('RMSE DFM')

subplot(2,1,2)
plot(xiRmseVAR, fRmseVAR, Color=colorVAR, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
title('RMSE VAR')
