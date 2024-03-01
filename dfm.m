%% MF Bayesian DFM with unbalanced data

close all
clc
seed=0;  
rng(seed);   

%% to delete
data = xlsread('rawdata2.xls',1,'r3:v860');     
y = data;  % includes Q and M variables -> mixed frequency
[T,N] = size(y);  % sample size and number of variables
%%

% data transformations for model
yMean = nanmean(y)';  % compute stats to add moments back after computing factor
ySd = nanstd(y)';  % avoid nan when computing stats
[yStandard,~,~] = standardMissObsOutl(y);  % standardization accounting for missing obs and outliers
yModel = yStandard';

% model specs
AR = 2;  % number of lags for each component
S =  N + N*Q + AR*M;   % number of states
m2q = [(1/3);(2/3);1;(2/3);(1/3)];  % month-to-quarter conversion
L = length(m2q);  % length of m2q conversion (additional lags from mapping) 


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
priorScaleErrCovIdio = 0.1*(priorShapeErrCovIdio-1);


%% Monte carlo

% specs
totDrawsKeep = 800;  % draws to keep 
totDrawsBurn = 200;  % burn-in
totDraws = totDrawsKeep + totDrawsBurn; 

% params (starting values)
% autoregressive coeff of common component (phi1, phi2)
coeffCommon = 0.5*ones(AR,1);

% error cov of common component 
errCovCommon = 1;  % normalized to one 

% loadings of common component (lambda or gamma) 
loadCommon = [ones(Q,1); 0.1*ones(M,1)];  % (first for quarterly, then for monthly)

% autoregressive coeff of idiosyncratic component (psi1, psi2)
coeffIdio = [0.2*ones(Q*AR,1); 0.2*ones(M*AR,1)];

% error cov of idiosyncratic component (sigma)
errCovIdio = 0.1*ones(Q+M,1);

% preallocate
posterCoeffCommon = zeros(length(coeffCommon), totDrawsKeep);
posterLoadCommon = zeros(length(loadCommon), totDrawsKeep);
posterCoeffIdio = zeros(length(coeffIdio), totDrawsKeep);
posterErrCovIdio = zeros(length(errCovIdio), totDrawsKeep);


%% posterior

% gibbs sampler
wait = waitbar(0, 'Gibbs');
tic
for itr = 1:totDraws
    
    % display iter
    disp(itr)
     
    % draw the latent, given the parameters and the data
    % =====================================================================

    [latent, yFill] = drawLatent(yModel, T, S, N, ...
        coeffCommon, errCovCommon, loadCommon, coeffIdio, errCovIdio, ...
        idx, Q, M, AR, m2q, L);
    latentCommon = latent(:,1);  % get latent of common component
    latentIdio = latent(:, [L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M]);  % get latent of idio component
    yEstimate = yFill * ySd(1) + yMean(1);  % add mean and std of target


    % draw the parameters, given the latent and the data
    % =====================================================================

    % common component
    %----------------------------------------------------------------------
    % draw the autoregressive coefficients 
    [coeffCommon] = drawCoeffCommon(latentCommon, priorCovCoeffCommon, ...
        priorMeanCoeffCommon, AR);

    % draw the factor loadings
    loadCommon(1,:) = drawLoadCommonMiss(Q, yModel, latentCommon, ...
        errCovIdio(Q), priorCovLoadCommon, priorMeanLoadCommon, ...
        idx, m2q, L);  % lower frequency   
    
    for jj = 1:M  % higher frequency
        loadCommon(Q+jj) = drawLoadCommon(Q+jj, yModel, latentCommon, ...
            coeffIdio(AR*Q+AR*(jj-1)+1:AR*Q+AR*jj), ...
            errCovIdio(Q+jj), ...
            priorCovLoadCommon, priorMeanLoadCommon, ...
            idx, AR);
    end        

    % idiosyncratic components
    %----------------------------------------------------------------------
    % draw the autoregressive coeff and error cov (both in same step)
    for jj = 1:Q  % lower frequency
        [errCovIdio(jj)] = drawErrCovIdio(jj, latentIdio, ...
            priorShapeErrCovIdio, priorScaleErrCovIdio, ...
            idx);
    end
    
    for jj = 1:M  % higher frequency
        [coeffIdio(AR*Q+AR*(jj-1)+1:AR*Q+AR*jj), errCovIdio(Q+jj)] = drawCoeffErrCovIdio(Q+jj, ...
            latentIdio, errCovIdio(Q+jj), priorCovCoeffIdio, ...
            priorMeanCoeffIdio, priorShapeErrCovIdio, priorScaleErrCovIdio, ... 
            idx, AR);
    end
 

    % collect draws for posterior empirical distributions
    % =====================================================================

    if itr > totDrawsBurn
        % autoregressive coeff in common component (phi1, phi2)
        posterCoeffCommon(:, itr-totDrawsBurn) = coeffCommon;
        % loadings in common components (lambda or gamma)
        posterLoadCommon(:, itr-totDrawsBurn) = loadCommon;
        % autoregressive coeff of idiosyncratic component (psi1, psi2)
        posterCoeffIdio(:, itr-totDrawsBurn) = coeffIdio;
        % error cov of idiosyncratic component (sigma)
        posterErrCovIdio(:, itr-totDrawsBurn) = errCovIdio;

        % factor (just common component not all latent)
        factorFiltered(:, itr-totDrawsBurn) = latentCommon;   
        %varfact(:,itr-totDrawsBurn) = errCovCommon;

        % predicted target (filled using transition equation in ss)
        yPredicted(:, itr-totDrawsBurn) = yEstimate;
    end

    waitbar(itr/totDraws, wait, sprintf('Gibbs sampling: %d%%', round(itr/totDraws*100)));

end
close(wait)
toc

% get summary of posteriors
factor = median(factorFiltered,2);
prct = [5 16 84 95];
factor_bands = prctile(factorFiltered, prct, 2);

target = median(yPredicted,2);


%% charts

% plot target along factor
time = 1952+2/12:1/12:2023+3/12;

figure;
subplot(2,1,1);
plot(time,target');
ylim([-15 10]);
xlim([time(1) time(end)]);
title('Monthly Real Activity Variable');

subplot(2,1,2);
plot(time,factor','-k');
hold on
plot(time,factor_bands(:,2:3)',':k');
xlim([time(1) time(end)]);
ylim([-15 10]);
title('Monthly Factor of Economic Activity');

figure;
plot(time,factor','r', LineWidth=1.5);
hold on
plot(time,target', 'b', LineWidth=1.5);
hold off
xlim([time(1) time(end)]);
ylim([-15 10]);
legend('factor', 'variables')
title('Monthly Observed vs Index of Economic Activity');

disp(corr(target,factor))

