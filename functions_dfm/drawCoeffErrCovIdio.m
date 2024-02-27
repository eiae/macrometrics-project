function [B,sigma2] = drawCoeffErrCovIdio(ndx,zz,sigma2,Sigma0,B0,T0,D0,index,AR)
%% Generate draw of autoregressive coeff and error cov in idiosyncratic component for monthly
% Input:
% - sample of states
% - prior hyperparams 
% - loadings in common component
% - error cov idiosyn
% Output:
% - draw of autoregressive coeff in idiosyncratic component
% -------------------------------------------------------------------------
     
% use selector vector for idiosyncratic components
z = zz(:,ndx);
Y = delIfCond(z,1-index(:,ndx));
X = zeros(length(Y),AR);

% AR estimation for coefficients (Bayesian OLS regression)
for ii = 1:AR
    X(:,ii) = lagMaker(Y,ii);
end

% remove missing obs
Y = Y(AR+1:end);
X = X(AR+1:end,:);

% additional terms for moments
TT = size(X,1);
     
% posterior moments for coeff|erro cov => Gaussian(posterior_mean, posterior_variance)
M = inv(Sigma0+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y);
V = inv(Sigma0+(1/sigma2)*(X'*X));

chck=-1;
while chck<0  % check for stability
    % posterior draw for coeff|error cov => Gaussian(posterior_mean, posterior_variance) 
    B = M + (randn(1,AR)*chol(V))';
    b = [B'; [eye(AR-1), zeros(AR-1,1)]];
    ee = max(abs(eig(b)));
    if ee<=1
        chck=1;
    end
end

% posterior moments for error cov|coeff => invGamma(posterior_shape, posterior_scale)
resids = Y - X*B;  % residuals
T1 = T0+TT/2;
D1 = (D0+(resids)'*(resids)/2);  % SSR

% posterior draw for error cov|coeff => invGamma(posterior_shape, posterior_scale)
sigma2 = 1/gamrnd(T1,1/D1);

end

