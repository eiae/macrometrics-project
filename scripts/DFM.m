%% MF Bayesian DFM with unbalanced data

% Disclaimer: code adapted from class which is based on routines written
% by Gabriel Pérez-Quirós and Danilo Leiva-León

% clear
clc

% data transformations for model
yMean = nanmean(yDFM)';  % compute stats to add moments back after computing factor
ySd = nanstd(yDFM)';  % avoid nan when computing stats
%[yStandard,~,~] = standardMissObsOutl(yDFM);  % standardization accounting for missing obs and outliers
[yStandard,~,~] = standardMissObs(yDFM);  % standardization accounting for missing obs and outliers
yModel = yStandard';

% model specs
AR = 2;  % number of lags for each component
m2q = [(1/3);(2/3);1;(2/3);(1/3)];  % month-to-quarter conversion
L = length(m2q);  % length of m2q conversion (additional lags from mapping) 
S =  L + L*Q + AR*M;   % number of states


%% priors

% autoregressive coeff of common component (phi1, phi2)
priorMeanCoeffCommon = zeros(AR,1);              
priorCovCoeffCommon = eye(AR);

% loadings of common component (lambda or gamma)
priorMeanLoadCommon = 0; 
priorCovLoadCommon = 1;

% autoregressive coeff of idiosyncratic components (psi1, psi2)
priorMeanCoeffIdio = zeros(AR,1); 
priorCovCoeffIdio = eye(AR);

% error cov of idiosyncratic components (sigma)
priorShapeErrCovIdio = 0;
priorScaleErrCovIdio = 0.1;  %0.1*(priorShapeErrCovIdio-1) -> negative which can create issues since posterior scale becomes negative and cannot draw from Gamma distribution


%% init simulation

% specs
%totDraws = totDrawsKeep + totDrawsBurn; 

% params (starting values)
% autoregressive coeff of common component (phi1, phi2)
coeffCommon = 0.5*ones(AR,1);

% error cov of common component 
errCovCommon = 1;  % normalized to one 

% loadings of idiosyncratic component (lambda or gamma) 
loadIdio = [ones(Q,1); 0.1*ones(M,1)];  % (first for quarterly, then for monthly)

% autoregressive coeff of idiosyncratic component (psi1, psi2)
coeffIdio = [0.2*ones(Q*AR,1); 0.2*ones(M*AR,1)];

% error cov of idiosyncratic component (sigma)
errCovIdio = 0.1*ones(Q+M,1);

% preallocate
posterCoeffCommon = zeros(length(coeffCommon), totDrawsKeep);
posterLoadCommon = zeros(length(loadIdio), totDrawsKeep);
posterCoeffIdio = zeros(length(coeffIdio), totDrawsKeep);
posterErrCovIdio = zeros(length(errCovIdio), totDrawsKeep);

factorFiltered = [];
yPredicted = [];
yPredictedCommon = [];

correlTargetPredAndPredComm = [];
correlTargetObsAndPredComm = [];

forecastTarget = [];


%% estimation (MCMC -> posterior)

% gibbs sampler
disp('Starting Mixed Frequency Bayesian DFM estimation...')
wait = waitbar(0);

%tic
for itr = 1:totDraws
    
    % display iter
    %disp(itr)
     
    % draw the latent, given the parameters and the data
    % =====================================================================

    [latent, yFill, yFillCommon, ...
        measureCoeff, measureCov, ...
        transCoeff, transCov] = drawLatent(yModel, T, S, N, ...
        coeffCommon, errCovCommon, loadIdio, coeffIdio, errCovIdio, ...
        idx, Q, M, AR, m2q, L);
    latentCommon = latent(:,1);  % get latent of common component
    latentIdio = latent(:, [L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M]);  % get latent of idio component
    %yEstimateNoMoments = yFill(:,1);
    yEstimate = yFill(:,1) * ySd(1) + yMean(1);  % add mean and std of target
    yEstimateCommon = yFillCommon(:,1) * ySd(1) + yMean(1);  % add mean and std of target

    % draw the parameters, given the latent and the data
    % =====================================================================

    % common component
    %----------------------------------------------------------------------
    % draw the autoregressive coefficients 
    coeffCommon = drawCoeffCommon(latentCommon, priorCovCoeffCommon, ...
        priorMeanCoeffCommon, AR);

    % idiosyncratic components
    %----------------------------------------------------------------------
    % draw the factor loadings
    for jj = 1:Q  % lower frequency   
        loadIdio(jj) = drawLoadIdioLF(jj, yModel, latentCommon, ...
            errCovIdio(jj), priorCovLoadCommon, priorMeanLoadCommon, ...
            idx, m2q, L);  
    end

    for jj = 1:M  % higher frequency
        loadIdio(Q+jj) = drawLoadIdioHF(Q+jj, yModel, latentCommon, ...
            coeffIdio(AR*Q+AR*(jj-1)+1:AR*Q+AR*jj), ...
            errCovIdio(Q+jj), ...
            priorCovLoadCommon, priorMeanLoadCommon, ...
            idx, AR);
    end        

    % draw the autoregressive coeff and error cov (both in same step)
    for jj = 1:Q  % lower frequency
        errCovIdio(jj) = drawErrCovIdioLF(jj, latentIdio, ...
            priorShapeErrCovIdio, priorScaleErrCovIdio, ...
            idx);
    end
    
    for jj = 1:M  % higher frequency
        [coeffIdio(AR*Q+AR*(jj-1)+1:AR*Q+AR*jj), errCovIdio(Q+jj)] = drawCoeffErrCovIdioHF(Q+jj, ...
            latentIdio, errCovIdio(Q+jj), priorCovCoeffIdio, ...
            priorMeanCoeffIdio, priorShapeErrCovIdio, priorScaleErrCovIdio, ... 
            idx, AR);
    end
 

    % collect draws for posterior empirical distributions
    % =====================================================================

    if itr > totDrawsBurn
        % autoregressive coeff in common component (phi1, phi2)
        posterCoeffCommon(:, itr-totDrawsBurn) = coeffCommon;
        % loadings in idiosyncratic components (lambda or gamma)
        posterLoadCommon(:, itr-totDrawsBurn) = loadIdio;
        % autoregressive coeff of idiosyncratic component (psi1, psi2)
        posterCoeffIdio(:, itr-totDrawsBurn) = coeffIdio;
        % error cov of idiosyncratic component (sigma)
        posterErrCovIdio(:, itr-totDrawsBurn) = errCovIdio;

        % factor (just common component not all latent)
        factorFiltered = [factorFiltered, latentCommon];   
        %varfact(:,itr-totDrawsBurn) = errCovCommon;

        % predicted target (filled using transition equation in ss)
        yPredicted = [yPredicted, yEstimate];

        % predicted target (filled using common part of transition equation in ss)
        yPredictedCommon = [yPredictedCommon, yEstimateCommon];

        % correlation
        correlTargetPredAndPredComm = [correlTargetPredAndPredComm; corr(yEstimate,yEstimateCommon)];
        
        nonMissTarget = ~isnan(yDFM(:,1));  
        yTargetObs = yDFM(nonMissTarget,1);
        yEstimateCommonQ = yEstimateCommon(nonMissTarget);
        correlTargetObsAndPredComm = [correlTargetObsAndPredComm; corr(yTargetObs,yEstimateCommonQ)];

        % forecast
        
        %yForecast(1:AR) = yEstimate(end-AR+1:end); % starting values for forecast based on AR lags 
        %B = [1; 0.3; 0.1];  % random walk process with 2 lags
        %cfactor = 1;  % sqrt(sigma2);  % standard deviation of the error cov
        % in the loop:
        % yForecast(h) = [1 yForecast(h-1) yForecast(h-2)]*B + (randn(1,1)*cfactor);
                
        yForecast = zeros(H+1,1);  % create new forecast each iteration
        yForecast(1) = yEstimate(end);
        latentForecast = zeros(H+1,S);
        latentForecast(1,:) = latent(end,:);
        stdError = sqrt(diag(transCov));  % standard deviation of the error cov in transition equation
        
        for h = 2:H+1  
            % transition equation to simulate forecast using VAR(1) dynamics
            latentForecast(h,:) = transCoeff^h * latentForecast(h-1,:)' + randn(S,S)*stdError; %.*stdError;  %randn(S,S)*stdError
            % measurement equation to map latent and observables
            yForecast(h) =  latentForecast(h,:) * measureCoeff(1,:)';  %measureCoeff(1,1:L)'
        end
        yForecastTarget = yForecast(:,1);
        forecastTarget = [forecastTarget, [yEstimate; yForecastTarget(2:end)]];  % concatenate history with forecast (taking into account initial values taken last observation depending on number of lags)
    end

    waitbar(itr/totDraws, wait, sprintf('Gibbs sampling for DFM: %d%%', round(itr/totDraws*100)));

end
close(wait)
%toc


%% results

% get observable target
obs = yDFM(:,1);

% get summary of posteriors
target = median(yPredicted,2);  % predicted with common and idio
common = median(yPredictedCommon,2);  % predicted with common 

factor = median(factorFiltered,2);
factor_bands = prctile(factorFiltered, prct, 2);

% correlation monthly predicted with idio + common vs only common
correlPred = median(correlTargetPredAndPredComm,1);
fprintf('Correlation between predicted target using common component vs common and idiosyncratic component (monthly):\n %4.2f \n', correlPred)

% correlation quarterly observable vs predicted common
correlObs = median(correlTargetObsAndPredComm,1);
fprintf('Correlation between target vs predicted target using common component (quarterly):\n %4.2f \n', correlObs)

% forecast
forecast = median(forecastTarget,2);
forecast_bands = prctile(forecastTarget, prct, 2);
datesForecast = [dates; (dates(end) + calmonths(1:H))'];

forecastPlot = [NaN(T-1,1); forecast(end-H:end)];
forecastPlotBands = [NaN(T-1,length(prct)); forecast_bands(end-H:end,:)];
historyPlot = [target; NaN(H,1)];


% collect quarterly predicted for comparison with target
datesQ = dates(~isnan(yDFM(:,1)));  % quarterly dates
targetQ = [];
commonQ = [];
obsQ = [];

i = 3*2;  % 3 periods per quarter and one lag due to growth
while i <= T  
    targetQ = [targetQ; target(i)];  
    commonQ = [commonQ; common(i)];  
    obsQ = [obsQ; obs(i)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end

% collect quarterly forecast for comparison to ECB
forecastQ = [];
forecast_bandsQ = [];

while i <= T+H  
    forecastQ = [forecastQ; forecast(i)];   
    forecast_bandsQ = [forecast_bandsQ; forecast_bands(i,:)];
    i = i+3;  % take only 1 in 3 values (3 month = 1 quarter)
end


% compute recursive mean to check MCMC convergence

% autoregressive coeff in common component (phi1, phi2)
recurMeanCoeffCommon = zeros(size(posterCoeffCommon,1), totDrawsKeep);
for j = 1:totDrawsKeep
    recurMeanCoeffCommon(:,j) = mean(posterCoeffCommon(:,1:j),2);
end

% loadings in idiosyncratic components (lambda or gamma)
recurMeanLoadCommon = zeros(size(posterLoadCommon,1), totDrawsKeep);
for j = 1:totDrawsKeep
    recurMeanLoadCommon(:,j) = mean(posterLoadCommon(:,1:j),2);
end

% autoregressive coeff of idiosyncratic component (psi1, psi2)
recurMeanCoeffIdio = zeros(size(posterCoeffIdio,1), totDrawsKeep);
for j = 1:totDrawsKeep
    recurMeanCoeffIdio(:,j) = mean(posterCoeffIdio(:,1:j),2);
end

% error cov of idiosyncratic component (sigma)
recurMeanErrCovIdio = zeros(size(posterErrCovIdio,1), totDrawsKeep);
for j = 1:totDrawsKeep
    recurMeanErrCovIdio(:,j) = mean(posterErrCovIdio(:,1:j),2);
end
