%% MF Bayesian DFM with unbalanced data

% specs
seed=0;  rng(seed);   % fix the seed if desired
idx = 1-isnan(y);  % matrix to select filled values

% transformations
yy_m = nanmean(y)';  % compute stats to add moments back after computing factor
yy_sd = nanstd(y)';  % avoid nan when computing stats

[Yb,~,~] = standardMissObsOutl(y);  % standardization accounting for missing obs and outliers
Y = Yb';

% model specs
AR = 2;  % number of lags for each component
S =  N + N*Q + AR*M;   % number of states
m2q = [(1/3);(2/3);1;(2/3);(1/3)];  % month-to-quarter conversion
L = length(m2q);  % length of m2q conversion lags for selector vectors 


%% priors

% for autoregressive coeff in common component (phi1, phi2)
T0 = [0;0];              
R0 = eye(2);

% for factor loadings in common components (lambda or gamma)
T00 = 0; 
R00 = 1;

% for autoregressive coeff in idiosyncratic component (psi1, psi2)
T0p = zeros(AR,1); 
R0p = eye(AR);

% for error cov of idiosyncratic component (sigma)
V0 = 0;
D0 = 0.1*(V0-1);


%% Monte carlo

% specs
N1 = 800;  % draws to keep 
N0 = 200;  % burn-in
Totiter = N0 + N1; 

% starting values
% for autoregressive coeff in common component (phi1, phi2)
PHI1TT = 0.6;  
PHI2TT = 0.1;

% for error cov of common component 
SIG2FTT = 1;  % normalized to one 

% for factor loadings in common components (lambda or gamma) 
LAM_TT = [ones(Q,1); ones(M,1)*.1];  % (1 for GDP, .1 for all else)

% for autoregressive coeff in idiosyncratic component (psi1, psi2)
PSI_TT = [0.2*ones(Q*AR,1); 0.2*ones(M*AR,1)];

% for error cov of idiosyncratic component (sigma)
SIG2_TT = 0.1*ones(Q+M,1);

% preallocate
lamda = zeros(length(LAM_TT), Totiter-N0);


%% posterior

% gibbs sampler
wait = waitbar(0, 'Gibbs');
tic
for itr=1:Totiter
    
    % display iter
    disp(itr)
     
    % 1) draw the states, given the parameters and the data
    [ZTTALL,Yfill] = drawLatent(Y,T,S,N,PHI1TT,PHI2TT,SIG2FTT,LAM_TT,PSI_TT,SIG2_TT,idx,Q,M,AR,m2q,L);
    FTT = ZTTALL(:,1);  % get common factor
    ITT = ZTTALL(:,[L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M]);  % individual selection in vector
    YHAT = Yfill*yy_sd(1)+yy_m(1);  % add mean and std of target


    % 2) draw the parameters, given the factor and the data
    % common component
    % draw the autoregressive coefficients 
    [PHI1TT,PHI2TT] = drawCoeffCommon(FTT,R0,T0,AR);

    % draw the factor loadings
    LAM_TT(1,:) = drawLoadCommonMiss(1,Y,FTT,SIG2_TT(1),R00,T00,idx,m2q,L);  % lower frequency   
    
    for iiv = 1:M  % higher frequency
        LAM_TT(Q+iiv) = ...
            drawLoadCommon(Q+iiv,Y,FTT,PSI_TT(AR*Q+AR*(iiv-1)+1:AR*Q+AR*iiv),SIG2_TT(Q+iiv),R00,T00,idx, AR);
    end        


    % idiosyncratic components
    % draw the autoregressive coeff and error cov (both in same step)
    for iiv = 1:Q  % lower frequency
        [SIG2_TT(iiv)] = drawErrCovIdio(iiv,ITT,V0,D0,idx);
    end
    
    for iiv = 1:M  % higher frequency
        [PSI_TT(AR*Q+AR*(iiv-1)+1:AR*Q+AR*iiv),SIG2_TT(Q+iiv)] = ...
            drawCoeffErrCovIdio(Q+iiv,ITT,SIG2_TT(Q+iiv),R0p,T0p,V0,D0,idx, AR);
    end
 

    % collect draws into empirical posterior distributions
    if itr>N0
        % for factor loadings in common components (lambda or gamma)
        lamda(:,itr-N0) = LAM_TT;
        
        % for factor (just common component not all states)
        fact(:,itr-N0) = FTT;   
        %varfact(:,itr-N0) = SIG2FTT;

        % for filled data
        forecastY(:,itr-N0) = YHAT;
    end

    waitbar(itr/Totiter, wait, sprintf('Gibbs sampling: %d%%', round(itr/Totiter*100)));

end
close(wait)
toc

% objects of interest
factor = median(fact,2);
factor_bands = prctile(fact,[5 16 84 95],2);
loadings = median(lamda')';
target = median(forecastY')';


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


