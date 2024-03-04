function [Rs,Qs,Hs,Fs] = stateSpaceMat(par, m2q, N, Q, M, S, L, AR)
%% Build state-space param matrices
% Input:
% - model params
% - variance of factor
% - month2quarter conversion vector
% - number of variables
% Output:
% - state-space coeff and innovation/error covariance matrices
% -------------------------------------------------------------------------

% automatic construction and preallocation of matrices

% H: coeff matrix obs equation
h1 = [m2q'*par(Q) m2q' zeros(1,AR*M)];  % lambda/gamma with m2q map (loading); AR*M since enlarge state vector according to AR (since state eq is VAR(1) in companion form)
h2 = alternateMat(M,AR*M,AR);
h3 = [par(Q+1:Q+M) zeros(M,L+(L-1)) h2];  % lambda/gamma (loading); L+(L-1) since you define first column of H' and then have 2*L (m2q for common and m2q for idio)
Hs = [h1; h3];      

% R: innovation covariance matrix obs equation
Rs = zeros(N,N);                   

% F: coeff matrix state equation
z1 = par(Q+M+1:Q+M+AR);  % phi (autoregressive of common)
z1 = z1*0;  % set to zero to avoid having the M2Q lags of quarterly idio
z2 = par(Q+M+AR+1:Q+M+AR*2);  % psi of quarterly (autoregressive of idio)
pamz3 = par(Q+M+AR*2+1:Q+M+AR*2+M*AR); % psi of monthly (autoregressive of idio)
z3 = alternateIdioMat(AR*M, AR*M, AR, pamz3);

f1 = [z1' zeros(1,S-AR); eye(L-1) zeros(L-1,S-(L-1))];                  
f2 = [zeros(1,L) z2' zeros(1,S-AR-L); zeros(L-1,L) eye(L-1) zeros(L-1,S-AR*2-L)];
f3 = [zeros(AR*M,L*2) z3];
Fs = [f1; f2; f3];

% Q: innovation covariance matrix state equation
parmq1 = par(Q+M+AR*2+M*AR+1:end);  % error variance of common and idio components
q1 = zeros(L,1);
q1(1) = parmq1(end);
q2 = zeros(L,1);
q2(1) = parmq1(1);
q3 = zeros(M*AR,1);
counter = 1;
for j = 1:M*AR
    if mod(j,2) == 1
        q3(j) = parmq1(counter+1);
        counter = counter + 1;
    end
end
q4 = [q1' q2' q3'];
% q4 = (q4).^2;  % no squaring since modeling variances (inv Gamma priors); if uncomment run into singularity -> solved with function <invpd()>?
Qs = diag(q4);
   
end
