function [Zdraw,yfilled] = drawLatent(vars,T,S,N,phi1,phi2,sig2fact,lam,psi,sig2,index,Q,M,AR,m2q,L)
%% Construct Kalman filter and Kalman smoother (KC algorithm)
% Input:
% - model params 
% Output:
% - sample/draws of states
% -------------------------------------------------------------------------

% specs
data = vars';  % data
params = [lam; phi1; phi2; psi; sig2; sig2fact];  % stack params (follow this order!)

% build state-space param matrices
[RR,QQ,H,F] = stateSpaceMat(params, m2q, N, Q, M, S, L, AR); 

% preallocate
Zp = zeros(S,1);      
Pp = eye(S);
Z = zeros(T,S);
P = zeros(S,S,T);


%% Kalman forward recursion (KF)

% filtering
for ii = 1:T
    % deal with missing values by reducing dimensions of state-space (Durbin-Koopman method)
    y = data(ii, index(ii,:)==1)';  % only use non-missing data for recursion 
    Hit = H(index(ii,:)==1,:);
    Rit = RR(index(ii,:)==1, index(ii,:)==1);
    
    % prediction error
    nu = y - Hit*Zp;  
    OM = Hit*Pp*Hit'+ Rit;    

    % Kalman gain
    invOM = inv(OM);  % invert separately for performance
    K = Pp*Hit'*invOM;
    
    % update step
    Ztt = Zp + K*nu;
    Ptt = Pp - K*Hit*Pp;
    
    % prediction step is different (smaller dimension)
    if ii < T
        Zp = F*Ztt;
        Pp = F*Ptt*F' + QQ;
    end
    
    % save iteration for smoother
    Z(ii,:) = Ztt';
    P(:,:,ii) = Ptt;
    
end
 
% this line here is again to match the previous ordering
% if AR==2 && M==4 && Q==1, yhat_corr=Z(:,[1,7:10])*H(:,[1,7:10])'; end


%% Kalman backward recursion (KC)

% reconstruct distribution of states to get first draw 
Zdraw(T,:) = mvnrndAlt(Ztt,Ptt,1);  % Zdraw(T|T) ~ N(S(T|T),P(T|T)); first iteration in smoother = last iteration in filter


% selection of matrices due to sparsity in the error covariance matrix in 
% state equation (efficiency gain and avoid singular matrix, see matrices5variables.pdf)
km = [1, L+1:L:L+L*Q, L+L*Q+1:AR:L+L*Q+AR*M];  % size of jumps in rows of matrix based on m2q conversion

% ordering
%if AR==2 && M==4 && Q==1, km = [1,2,3,4,5,6]; end

Qstar = QQ(km,km);
Fstar = F(km,:);

% smoothing
for ii = T-1:-1:1
    
    % define state mean and variance for backward recursion
    Zf = Zdraw(ii+1,km)';  % ht+1
    Ztt = Z(ii,:)';  % ht|t (= ht|T)
    Ptt = P(:,:,ii);  % Pt|t (= Pt|T)
    
    % update step
    % define variables that simplify computation of moments
    OMsmooth = Fstar*Ptt*Fstar' + Qstar;  % prediction variance state eq
    Ksmooth = Ptt*Fstar'*inv(OMsmooth);  % Kalman gain
    nusmooth = Zf - Fstar*Ztt;  % state error

    % state moments to draw state
    Zmean = Ztt + Ksmooth*nusmooth;  % ht+1|t, state mean (see Bayesian.pdf or eq. 13.6.11 and 13.6.12 Hamilton)
    Zvar = Ptt - Ksmooth*Fstar*Ptt;  % Pt+1|t, state variance    
    Zdraw(ii,:) = mvnrndAlt(Zmean,Zvar,1);  % draw

end
    
% ordering
% if AR==2 && M==4 && Q==1, Fdraw=Zdraw(:,d_to_e); end

%% fill missing obs
% use transition equation for forecasting missing obs or matching data
yfilled = Z*H(1,:)';  
     
end
        
