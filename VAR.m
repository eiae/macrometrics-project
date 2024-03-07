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
target = yVAR(:,1);
pred = median(yPredicted, 3); 


%% charts 

% levels
figure;
plot(dates, target, 'Linewidth',2, 'marker','o', 'color','r'); 
hold on; 
plot(dates, pred, 'Linewidth',2, 'color', 'b')
legend('observed (quarterly)', 'index (monthly)')
axis tight
title('Quarterly Level of Observed vs Index of Economic Activity');

% growth rates
targetGrowth = 100*( target(L+1:end,1) - target(1:end-L,1) );
predGrowth = 100*( pred(L+1:end,1) - pred(1:end-L,1) );

figure;
plot(dates(L+1:end), targetGrowth, 'Linewidth',2, 'marker','o', 'color','r');
hold on; 
plot(dates(L+1:end), predGrowth, 'Linewidth',2, 'color', 'b')
legend('observed (quarterly)', 'index (monthly)')
axis tight
title('Quarterly Growth of Observed vs Index of Economic Activity');
