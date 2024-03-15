%% MF Bayesian VAR with unbalanced data

% Disclaimer: code adapted from Empirical Macro Toolbox
% by Filippo Ferroni and Fabio Canova
% https://github.com/naffe15/BVAR_/tree/master

close all
clc
seed=0;  
rng(seed);  
 
% model specs
lags = 2;  % number of lags
options.mf_varindex = 1;  % position of target variable
options.K = 1000;  % total iterations
options.priors.name = 'Minnesota';  % prior type with standard hyperparams
options.fhor = H; 
options.noprint = 1;

% given we are working with level data the M2Q match is
% Yq(t) = 1/3*Ym(t) + 1/3*Ym(t-1) + 1/3*Ym(t-2)
m2q = [(1/3);(1/3);(1/3)];  % month-to-quarter conversion
L = length(m2q);  % length of m2q conversion (additional lags from mapping) 


%% estimation

% estimate BVAR
bvarmf = bvar_(yVAR, lags, options);  

% sort predicted endog vector
yPredicted = sort(bvarmf.yfill, 3);
yPredicted = yPredicted(:,1,:);  % select target variable 

% get summary of posteriors
obsAlt = yVAR(:,1);
targetAlt = median(yPredicted, 3); 

% forecast 
forecastTargetAlt = squeeze(bvarmf.forecasts.with_shocks(:,1,:));

forecastAlt = median(forecastTargetAlt, 2);
forecast_bandsAlt = prctile(forecastTargetAlt, prct, 2);

datesForecast = [dates; (dates(end) + calmonths(1:H))'];

forecastPlotAlt = [NaN(T,1); forecastAlt];
forecastPlotBandsAlt = [NaN(T,length(prct)); forecast_bandsAlt];
historyPlotAlt = [targetAlt; NaN(H,1)];


%% charts 

% levels
figure;
plot(dates, obsAlt, 'Linewidth',2, 'marker','o', 'color','r'); 
hold on; 
plot(dates, targetAlt, 'Linewidth',2, 'color', 'b')
legend('observed (quarterly)', 'index (monthly)')
axis tight
title('Quarterly Level of Observed vs Index of Economic Activity');

% growth rates
obsGrowth = 100*( obsAlt(L+1:end,1) - obsAlt(1:end-L,1) );
targetGrowth = 100*( targetAlt(L+1:end,1) - targetAlt(1:end-L,1) );

figure;
plot(dates(L+1:end), obsGrowth, 'Linewidth',2, 'marker','o', 'color','r');
hold on; 
plot(dates(L+1:end), targetGrowth, 'Linewidth',2, 'color', 'b')
legend('observed (quarterly)', 'index (monthly)')
axis tight
title('Quarterly Growth of Observed vs Index of Economic Activity');


% plot point and density forecast
figure;
plot(datesForecast, historyPlotAlt, Color='b', Linewidth=2)
hold on;
plot(datesForecast, forecastPlotAlt, Color='r', Linewidth=1.5, LineStyle='-')
hold on;
plot(datesForecast(T:end)' , [forecastPlotBandsAlt(T:end,1), forecastPlotBandsAlt(T:end,end)],...
    Color='r', LineWidth=1, LineStyle=':')
hFill = [datesForecast(T:end)' fliplr(datesForecast(T:end)')];
inBetween = [forecastPlotBandsAlt(T:end,1)', fliplr(forecastPlotBandsAlt(T:end,end)')];
fill(hFill , inBetween, 'r', FaceAlpha=0.2, LineStyle='none');
legend('history', 'point forecast', sprintf('credible bands %d%%',prct(end)))
axis tight
grid on
title('Forecast of Monthly Predicted Real Activity Variable');
