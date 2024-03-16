%% Comparison analysis between models

% Disclaimer: code written
% by Erik Andres Escayola

% clear
clc



%%

% redefining some variables for the exercise
H = 12;
Tfull = T; 
min_window_size = Tfull-12*2;
RmseDFM = zeros(Tfull-H, 1);

yDFMFull = preprocessDFM(yraw,selected,Tfull,N,Q);
obsFull = yDFMFull(:,1);
obsFullQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= T  
    obsFullQ = [obsFullQ; obsFull(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% progressively expanding window and reestimating model
counter = 1;
for i = min_window_size:Tfull-H

    % select subsample
    subsample = yraw(1:i, :);
    T = size(subsample,1);
    idx = 1-isnan(subsample);
    dates = dates(1:i,:);

    % transformations DFM -> growth rates
    yDFM = preprocessDFM(subsample,selected,T,N,Q);

    % run DFM
    run(fullfile(workpath,'DFM.m'))
    %run(fullfile(workpath,'DFM_CHART.m'))

    % compare forecast with actual data
    % CHECK PERIODS!!!!
    actual_value = obsFullQ(end-H/3+1:end);
    forecast_value = forecastQ; %forecastQ(end-H/3+1:end);
    
    % compute rmse
    %error = actual_value - forecast_value;
    %RMSE= sqrt(mean(error.^2))
    RMSE = rmse(actual_value, forecast_value);
    
    % store RMSE value
    RmseDFM(counter) = RMSE;

    counter = counter + 1;
end

% Plot RMSE values
plot(RmseDFM)
xlabel('Number of observations')
ylabel('Root Mean Squared Error (RMSE)')
title('RMSE DFM')