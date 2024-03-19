function [yForecast, latentForecast] = forecastSparseCov(yEstimate, latent, transCov, transCoeff, measureCoeff, Q, M, S, L, AR, H)
%% Forecast algorithm when transition error covariance is sparse
% Input:
% - filled target
% - filtered latent (draw)
% - transition error covariance
% - transition coeff matrix
% - measurement coeff matrix
% - number of quarterly, monthly variables, latent variables, length of m2q
%   mapping, lags, forecast horizon
% Output:
% - target forecast
% - latent variable forecast
% -------------------------------------------------------------------------

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
transCoeffRed = transCoeff(pick,:);  % select all columns to have all AR dynamics

% forecast by simulating the economy on the reduced transition
% equation (contains all dynamics), then recover original size of
% objects, finally match latent with observables in measurement
% equation
for h = 2:H+1  
    % reduced transition equation to simulate forecast using companion VAR(1) dynamics
    latentForecast(h,pick) = transCoeffRed * latentForecast(h-1,:)' + transStdRed*randn(sizePick,1);  % error/shock is a standard random normal with covariance matrix equal the transition error covariance

    % reconstruct latent vector with original size:
    % we have 5 cases due to the m2q mapping which contains the 
    % current month + four lags of the months to reconstruc the 
    % quarter. Yet, in each of the first 5 horizons there is a 
    % dynamic combination between previous forecasts and initial 
    % forecast based on last draw of the latent vector
    if h==2       
        % update subset of elements with initial forecast
        latentForecast(h,pickAlt) = latentForecastInit(pickAlt);
    
    elseif h==3        
        % update subset of elements with 1-horizon previous forecast
        pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];   %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev) = latentForecast(h-1,pickPrev-1); 
       
        % update subset of elements with initial forecast
        pickAlt = [L-2:L, L+L*Q-2:L+L*Q];  % here not automatized for Q>1  %[3,4,5,8,9,10] for 5 variables (1Q, 4M)
        latentForecast(h,pickAlt) = latentForecastInit(pickAlt-2);

    elseif h==4      
        % update subset of elements with 1-horizon previous forecast
        pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  
        latentForecast(h,pickPrev) = latentForecast(h-1,pickPrev-1);
        
        % update subset of elements with 2-horizons previous forecast
        pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev2) = latentForecast(h-2, pickPrev2-2);
        
        % update subset of elements with initial forecast
        pickAlt = [L-1:L, L+L*Q-1:L+L*Q];  % here not automatized for Q>1  %[4,5,9,10] for 5 variables (1Q, 4M)
        latentForecast(h,pickAlt) = latentForecastInit(pickAlt-3);

    elseif h==5
        % update subset of elements with 1-horizon previous forecast
        pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev) = latentForecast(h-1,pickPrev-1);
        
        % update subset of elements with 2-horizons previous forecast
        pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev2) = latentForecast(h-2,pickPrev2-2);
        
        % update subset of elements with 3-horizons previous forecast    
        pickPrev3 = [L-1, L+L*Q-1];  % here not automatized for Q>1  %[4,9] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev3) = latentForecast(h-3,pickPrev3-3);
        
        % update subset of elements with initial forecast
        pickAlt = [L, L+L*Q];  % here not automatized for Q>1  %[5,10]for 5 variables (1Q, 4M)
        latentForecast(h,pickAlt) = latentForecastInit(pickAlt-4);  

    else       
        % update subset of elements with 1-horizon previous forecast
        pickPrev = [L-3, L+L*Q-3, (L+L*Q+1+AR-1):AR:(L+L*Q+AR*M)];  %[2,7,12,14,16,18] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev) = latentForecast(h-1,pickPrev-1);
        
        % update subset of elements with 2-horizons previous forecast
        pickPrev2 = [L-2, L+L*Q-2];  % here not automatized for Q>1  %[3,8] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev2) = latentForecast(h-2,pickPrev2-2);
        
        % update subset of elements with 3-horizons previous forecast
        pickPrev3 = [L-1, L+L*Q-1];  % here not automatized for Q>1  %[4,9] for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev3) = latentForecast(h-3,pickPrev3-3);
        
        % update subset of elements with 4-horizons previous forecast
        pickPrev4 = [L, L+L*Q];  % here not automatized for Q>1  %[5,10]for 5 variables (1Q, 4M)
        latentForecast(h,pickPrev4) = latentForecast(h-4,pickPrev4-4);
    end
    
    % measurement equation to map latent and observables
    yForecast(h) =  latentForecast(h,:) * measureCoeff(1,:)';  % only match target variable  % measureCoeff(1,1:L)'

end



end