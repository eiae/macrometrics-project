%% MF Bayesian DFM with unbalanced data

% Disclaimer: code adapted from class which is based on routines written
% by Gabriel Pérez-Quirós and Danilo Leiva-León


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
coeffIdio = [0.2*ones(Q*AR,1)*0; 0.2*ones(M*AR,1)];  % multiply by zero since low frequency variable has no idiosyncratic component

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

        % preallocate
        yForecast = zeros(H+1,1);  % create new forecast each iteration
        yForecast(1) = yEstimate(end);  % due to companion VAR(1) dynamics need to initialize forecast with last estimate (yet will only keep t+1 -> t+H+1 elements)
        
        latentForecastInit = latent(end,:);  % save last draw of latent vector for transition equation forecast
        
        latentForecast = zeros(H+1,S);  % create new forecast each iteration
        latentForecast(1,:) = latentForecastInit;  % due to companion VAR(1) dynamics need to initialize forecast with last estimate (yet will only keep t+1 -> t+H+1 elements)

        % define reduced matrices for simulation (avoid sparse transition covariance - Kim-Nelson section 8.2)
        pick = [1, L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M];  %[1,6,11,13,15,17] for 5 variables (1Q, 4M)
        sizePick = length(pick);

        pickAll = 1:S;
        checkPick = ~ismember(pickAll, pick);
        pickAlt = pickAll(checkPick);  % complement of the picker
        
        % first reduce dimensionality of transition equation for simulation 
        % then recover original size of objects (dynamically)

        % compute standard error of the shock to simulate the economy
        transCovRed = transCov(pick,pick);  % reduced transition equation error covariance to avoid sparse matrix
        transStdRed = chol(transCovRed)';  % take cholesky as to get a measure of standard error

        % define reduced transition coeff matrix
        transCoeffRed = transCoeff(pick,pick);

        % define reduced latent vector
        latentForecastRed = zeros(H+1,sizePick);
        latentForecastRed(1,:) = latentForecastInit(pick);  % due to companion VAR(1) dynamics need to initialize forecast with last estimate (yet will only keep t+1 -> t+H+1 elements)

        % forecast by simulating the economy on the reduced transition
        % equation (contains all dynamics), then recover original size of
        % objects, finally match latent with observables in measurement
        % equation
        for h = 2:H+1  
            % reduced transition equation to simulate forecast using companion VAR(1) dynamics
            latentForecastRed(h,:) = transCoeffRed * latentForecastRed(h-1,:)' + transStdRed*randn(sizePick,1);  % error/shock is a standard random normal with covariance matrix equal the transition error covariance
            
            % reconstruct latent vector with original size:
            % we have 5 cases due to the m2q mapping which contains the 
            % current month + four lags of the months to reconstruc the 
            % quarter. Yet, in each of the first 5 horizons there is a 
            % dynamic combination between previous forecasts and initial 
            % forecast based on last draw of the latent vector
            if h==2
                % update subset of elements with current forecast
                latentForecast(h,pick) = latentForecastRed(h,:);
                
                % update subset of elements with initial forecast
                latentForecast(h,pickAlt) = latentForecastInit(pickAlt);
            
            elseif h==3
                % update subset of elements with current forecast
                latentForecast(h,pick) = latentForecastRed(h,:);
                
                % update subset of elements with 1-horizon previous forecast
                pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];   %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev) = latentForecastRed(h-1,:);  % note that the reduced has already the correct picked positions ([1,6,11,13,15,17])
               
                % update subset of elements with initial forecast
                pickAlt = [L-2:L, L+L*Q-2:L+L*Q];  % here not automatized for Q>1  %[3,4,5,8,9,10] for 5 variables (1Q, 4M)
                latentForecast(h,pickAlt) = latentForecastInit(pickAlt-2);

            elseif h==4
                % update subset of elements with current forecast
                latentForecast(h,pick) = latentForecastRed(h,:);        
                
                % update subset of elements with 1-horizon previous forecast
                pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  
                latentForecast(h,pickPrev) = latentForecastRed(h-1,:);
                
                % update subset of elements with 2-horizons previous forecast
                pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev2) = latentForecastRed(h-2, pickPrev2-2);
                
                % update subset of elements with initial forecast
                pickAlt = [L-1:L, L+L*Q-1:L+L*Q];  % here not automatized for Q>1  %[4,5,9,10] for 5 variables (1Q, 4M)
                latentForecast(h,pickAlt) = latentForecastInit(pickAlt-3);

            elseif h==5
                % update subset of elements with current forecast
                latentForecast(h,pick) = latentForecastRed(h,:);
                
                % update subset of elements with 1-horizon previous forecast
                pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev) = latentForecastRed(h-1,:);
                
                % update subset of elements with 2-horizons previous forecast
                pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev2) = latentForecastRed(h-2,pickPrev2-2);
                
                % update subset of elements with 3-horizons previous forecast    
                pickPrev3 = [L-1, L+L*Q-1];  % here not automatized for Q>1  %[4,9] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev3) = latentForecastRed(h-3,pickPrev3-3);
                
                % update subset of elements with initial forecast
                pickAlt = [L, L+L*Q];  % here not automatized for Q>1  %[5,10]for 5 variables (1Q, 4M)
                latentForecast(h,pickAlt) = latentForecastInit(pickAlt-4);  

            else
                % update subset of elements with current forecast
                latentForecast(h,pick) = latentForecastRed(h,:);
                
                % update subset of elements with 1-horizon previous forecast
                pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev) = latentForecastRed(h-1,:);
                
                % update subset of elements with 2-horizons previous forecast
                pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev2) = latentForecastRed(h-2,pickPrev2-2);
                
                % update subset of elements with 3-horizons previous forecast
                pickPrev3 = [L-1, L+L*Q-1];  % here not automatized for Q>1  %[4,9] for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev3) = latentForecastRed(h-3,pickPrev3-3);
                
                % update subset of elements with 4-horizons previous forecast
                pickPrev4 = [L, L+L*Q];  % here not automatized for Q>1  %[5,10]for 5 variables (1Q, 4M)
                latentForecast(h,pickPrev4) = latentForecastRed(h-4,pickPrev4-4);
            end
            
            % measurement equation to map latent and observables
            yForecast(h) =  latentForecast(h,:) * measureCoeff(1,:)';  % only match target variable  % measureCoeff(1,1:L)'
        
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

i = 3*2;  % 3 periods per quarter and one lag due to growth
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
