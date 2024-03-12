function B = drawLoadIdioHF(ndx,yy,zz,psi,sigma2,Sigma0,B0,index,AR)
%% Generate draw of loadings in common component for monthly
% Input:
% - sample of states
% - prior hyperparams 
% - autoregressive coeff idiosyn
% - error cov idiosyn
% Output:
% - draw of loadings in common component
% -------------------------------------------------------------------------

yy0 = yy(ndx,:)';  % take target variable

% avoid missing in linear regression
yy=delIfCond(yy0,1-index(:,ndx));
ztt=delIfCond(zz,1-index(:,ndx));

% remove missing obs
Y = yy(AR+1:end,:); 
X = ztt(AR+1:end,:);

% rewrite the model in terms of AR polynomial
for ii=1:AR
    Y = Y - psi(ii)*yy(AR+1-ii:end-ii,:);  % endog variables
    X = X - psi(ii)*ztt(AR+1-ii:end-ii,:);  % exog variables = factor
end

% posterior moments for coeff|erro cov => Gaussian(posterior_mean, posterior_variance)
M = inv(Sigma0 + (1/sigma2)*(X'*X))*(inv(Sigma0)*B0 + (1/sigma2)*X'*Y);
V = inv(Sigma0 + (1/sigma2)*(X'*X));

% posterior draw for coeff|error cov => Gaussian(posterior_mean, posterior_variance)
B = M+(randn(1,1)*chol(V))';
 
end
