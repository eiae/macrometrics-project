%% Forecast exercise for VAR

% Disclaimer: code written
% by Erik Andres Escayola


%% forecast comparison between VAR and actual

% preallocate
RmseVAR = zeros(forecastPeriods, 1);
RmsedatesAlt = [];

% preprocess data and get actual with full sample
yVARFull = preprocessVAR(yraw,selected,Tfull,N,Q);
obsAltFull = [NaN(LAlt,1); 100*( yVARFull(LAlt+1:end,1) - yVARFull(1:end-LAlt,1) )];  % quarter-on-quarter growth rates based on monthly series of quartely data

datesfullAltQ = datesfull(~isnan(obsAltFull)); 
obsAltFullQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= Tfull  
    obsAltFullQ = [obsAltFullQ; obsAltFull(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% estimate model in expanding window 
counter = 1;
for window = minWindow:3:Tfull-H
    
    fprintf('VAR forecast iteration progress: %0.1f%%\n',  round(counter/forecastPeriods*100,1))

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
        % given we are using MF with the target variable as regressor, 
    % it only makes sense to compute RMSE for the periods of the forecast 
    % horizons (not in-sample since the data will be matched)
    backShiftAlt = length(forecastAltQ(1:end-HQ));
    actualValueVAR = obsAltFullQ(backShiftAlt+1:backShiftAlt+HQ);
    forecastValueVAR = forecastAltQ(end-HQ+1:end);


    disp('Actual vs forecast values for horizons:')
    disp([actualValueVAR,forecastValueVAR])
    
    % compute rmse
    %error = actual_value - forecast_value;
    %RMSE= sqrt(mean(error.^2))
    RmseValueVAR = rmse(actualValueVAR, forecastValueVAR);
    
    % store rmse value
    RmseVAR(counter) = RmseValueVAR;
    RmsedatesAlt = [RmsedatesAlt; datesfullAltQ(backShiftAlt)];

    counter = counter + 1;
end

% compute kernel density
[fRmseVAR, xiRmseVAR] = ksdensity(RmseVAR);