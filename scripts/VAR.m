%% MF Bayesian VAR with unbalanced data

% Disclaimer: code adapted from Empirical Macro Toolbox
% by Filippo Ferroni and Fabio Canova
% https://github.com/naffe15/BVAR_/tree/master

% clear
clc
 
% model specs
lags = AR;  % number of lags
options.mf_varindex = 1;  % position of target variable
options.priors.name = 'Minnesota';  % prior type with standard hyperparams
options.fhor = H; % horizon to forecast
options.noprint = 1;
%totDraws = totDrawsKeep + totDrawsBurn; 
options.K = totDraws;  % total iterations


% given we are working with level data the M2Q match is
% Yq(t) = 1/3*Ym(t) + 1/3*Ym(t-1) + 1/3*Ym(t-2)
m2q = [(1/3);(1/3);(1/3)];  % month-to-quarter conversion
L = length(m2q);  % length of m2q conversion (additional lags from mapping) 


%% estimation (MCMC -> posterior)

% gibbs sampler
bvarmf = bvar_(yVAR, lags, options); 

% collect objects of interest
posterCoeff = bvarmf.Phi_draws(:,:,end-totDrawsKeep+1:end);  % autoregressive coeff (Phi)
posterErrCov = bvarmf.Sigma_draws(:,:,end-totDrawsKeep+1:end);  % error cov (Sigma)

yPredictedTemp = sort(bvarmf.yfill, 3);  % sort predicted endog vector
yPredictedLvl = squeeze(yPredictedTemp(:,1,end-totDrawsKeep+1:end));  % select target variable 

forecastTemp = squeeze(bvarmf.forecasts.with_shocks(:,1,end-totDrawsKeep+1:end));
forecastTargetLvl = [yPredictedLvl; forecastTemp];  % join history and forecast


%% results

% get observable target and transform to growth rates
obsAlt = [NaN(L,1); 100*( yVAR(L+1:end,1) - yVAR(1:end-L,1) )];  % quarter-on-quarter growth rates based on monthly series of quartely data
yPredictedAlt = [NaN(L,totDrawsKeep); 100*( yPredictedLvl(L+1:end,:) - yPredictedLvl(1:end-L,:) )];  % add NaN to match sizes (since lose obs equal to frequency-relation length, M2Q = 3)
forecastTargetAlt = [NaN(L,totDrawsKeep); 100*( forecastTargetLvl(L+1:end,:) - forecastTargetLvl(1:end-L,:) )];  

% get summary of posteriors
targetAlt = median(yPredictedAlt, 2); 

% forecast 
forecastAlt = median(forecastTargetAlt, 2);
forecast_bandsAlt = prctile(forecastTargetAlt, prct, 2);

datesForecast = [dates; (dates(end) + calmonths(1:H))'];

forecastPlotAlt = [NaN(T-1,1); forecastAlt(end-H:end)];
forecastPlotBandsAlt = [NaN(T-1,length(prct)); forecast_bandsAlt(end-H:end,:)];
historyPlotAlt = [targetAlt; NaN(H,1)];

% levels (recall BVAR enters in levels)
obsAltLvl = yVAR(:,1);
targetAltLvl = median(yPredictedLvl, 2); 

% collect quarterly predicted for comparison with target
targetAltQ = [];
obsAltQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= T
    targetAltQ = [targetAltQ; targetAlt(i)];
    obsAltQ = [obsAltQ; obsAlt(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end


% collect quarterly forecast for comparison to ECB
forecastAltQ = [];
forecast_bandsAltQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= T+H  
    forecastAltQ = [forecastAltQ; forecastAlt(i)];   
    forecast_bandsAltQ = [forecast_bandsAltQ; forecast_bandsAlt(i,:)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end


% compute recursive mean to check MCMC convergence

% autoregressive coeff (Phi)
recurMeanCoeff = zeros(size(posterCoeff,1), N, totDrawsKeep);
for i = 1:N
    for j = 1:totDrawsKeep
        recurMeanCoeff(:,i,j) = mean(posterCoeff(:,i,1:j),3);
    end
end

% error cov (Sigma)
recurMeanErrCov = zeros(size(posterErrCov,1), N, totDrawsKeep);
for i = 1:N
    for j = 1:totDrawsKeep
        recurMeanErrCov(:,i,j) = mean(posterErrCov(:,i,1:j),3);
    end
end
