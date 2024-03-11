function [latentDraw,yfilled,yfilledFactor] = drawLatent(vars,T,S,N,phi,sig2fact,lam,psi,sig2,index,Q,M,AR,m2q,L)
%% Construct Kalman filter and Kalman smoother (KC algorithm)
% Input:
% - model params 
% Output:
% - sample/draws of states
% - filled target using factor and idio
% - filled target using factor only (to see explanatory power of model)
% -------------------------------------------------------------------------

% specs
data = vars';  % data
params = [lam; phi; psi; sig2; sig2fact];  % stack params (follow this order!)

% build state-space param matrices
[RR,QQ,H,F] = stateSpaceMat(params, m2q, N, Q, M, S, L, AR); 

% preallocate
h = zeros(S,1);      
P = eye(S);
hSmooth = zeros(T,S);
PSmooth = zeros(S,S,T);


%% Kalman forward recursion (KF)

% filtering
for ii = 1:T
    % deal with missing values by reducing dimensions of state-space (Durbin-Koopman method)
    y = data(ii, index(ii,:)==1)';  % only use non-missing data for recursion 
    Hit = H(index(ii,:)==1,:);
    Rit = RR(index(ii,:)==1, index(ii,:)==1);
    
    % prediction error
    nu = y - Hit*h;  
    O = Hit*P*Hit'+ Rit;    

    % Kalman gain
    invO = invNonPdMat(O);  % invert separately for performance
    K = P*Hit'*invO;
    
    % update step
    htt = h + K*nu;
    Ptt = P - K*Hit*P;
    
    % prediction step is different (smaller dimension)
    if ii < T
        h = F*htt;
        P = F*Ptt*F' + QQ;
    end
    
    % save iteration for smoother
    hSmooth(ii,:) = htt';
    PSmooth(:,:,ii) = Ptt;
    
end
 

%% Kalman backward recursion (KC)

% reconstruct distribution of states to get first draw 
latentDraw(T,:) = mvnrndAlt(htt,Ptt,1);  % stateDraw(T|T) ~ N(state(T|T), stateUncertainty(T|T)); first iteration in smoother = last iteration in filter

% selection of matrices due to sparsity in the error covariance matrix in 
% state equation (efficiency gain and avoid singular matrix, see matrices5variables.pdf)
pick = [1, L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M];  % size of jumps in rows of matrix based on m2q conversion
Fstar = F(pick,:);
Qstar = QQ(pick,pick);

% smoothing
for ii = T-1:-1:1
    
    % define state mean and variance for backward recursion
    hLead = latentDraw(ii+1,pick)';  % ht+1
    httSmooth = hSmooth(ii,:)';  % ht|t (= ht|T)
    PttSmooth = PSmooth(:,:,ii);  % Pt|t (= Pt|T)
    
    % update step
    % define variables that simplify computation of moments
    OSmooth = Fstar*PttSmooth*Fstar' + Qstar;  % prediction variance state eq
    KSmooth = PttSmooth*Fstar'*invNonPdMat(OSmooth);  % Kalman gain
    nuSmooth = hLead - Fstar*httSmooth;  % state error

    % state moments to draw state
    latentMean = httSmooth + KSmooth*nuSmooth;  % ht+1|t, state mean (see Bayesian.pdf or eq. 13.6.11 and 13.6.12 Hamilton)
    latentCov = PttSmooth - KSmooth*Fstar*PttSmooth;  % Pt+1|t, state variance    
    latentDraw(ii,:) = mvnrndAlt(latentMean,latentCov,1);  % draw

end

%% fill missing obs
% use transition equation for forecasting missing obs or matching data
yfilled = hSmooth*H';  
yfilledFactor = hSmooth(:,1:L)*H(:,1:L)';
     
end
        
