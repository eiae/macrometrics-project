%% Forecast exercise for DFM

% Disclaimer: code written
% by Erik Andres Escayola


%% forecast comparison between DFM and actual

% preallocate
RmseDFM = zeros(forecastPeriods, 1);
Rmsedates = [];

% preprocess data and get actual with full sample
yDFMFull = preprocessDFM(yraw,selected,Tfull,N,Q);
obsFull = yDFMFull(:,1);

datesfullQ = datesfull(~isnan(obsFull)); 
obsFullQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= Tfull  
    obsFullQ = [obsFullQ; obsFull(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% estimate model in expanding window 
counter = 1;
for window = minWindow:3:Tfull-H
    
    fprintf('DFM forecast iteration progress: %0.1f%%\n',  round(counter/forecastPeriods*100,1))

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
    
    disp('Actual vs forecast values for horizons:')
    disp([actualValueDFM,forecastValueDFM])

    % compute rmse
    %error = actual_value - forecast_value;
    %RMSE= sqrt(mean(error.^2))
    RmseValueDFM = rmse(actualValueDFM, forecastValueDFM);

    % store rmse value
    RmseDFM(counter) = RmseValueDFM;
    Rmsedates = [Rmsedates; datesfullQ(backShift)];

    counter = counter + 1;
end

% compute kernel density
[fRmseDFM, xiRmseDFM] = ksdensity(RmseDFM);

