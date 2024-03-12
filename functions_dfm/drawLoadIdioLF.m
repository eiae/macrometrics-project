function B = drawLoadIdioLF(ndx,yy,z,sigma2,Sigma0,B0,index0,m2q,L)
%% Generate draw of loadings in common component for quarterly 
% Input:
% - sample of states
% - prior hyperparams 
% - autoregressive coeff idiosyn
% - error cov idiosyn
% Output:
% - draw of loadings in common component
% -------------------------------------------------------------------------

yy = yy(ndx,:)';  % take target variable
y = yy(L:end,:);  % adjust size to account for m2q lags
X = m2q(1)*z(L:end) + m2q(2)*z((L-1):end-1) + m2q(3)*z((L-2):end-2) + ...
    m2q(4)*z((L-3):end-3) + m2q(5)*z((L-4):end-4);
index = index0(L:end,:);  % adjust size to account for m2q lags

% avoid missing in linear regression
y = delIfCond(y, 1-index(:,ndx));
X = delIfCond(X, 1-index(:,ndx));

% given the assumption of idiosyncractic component have "iid errors", 
% need to account M2Q relation when computing variance
Sigma2 = sigma2*(sum(m2q.^2));  % (1/3)^2 + (2/3)^2 + ... + (1/3)^2

% posterior moments for coeff|erro cov => Gaussian(posterior_mean, posterior_variance)
M = inv(Sigma0 + (1/Sigma2)*(X'*X))*(inv(Sigma0)*B0 + (1/Sigma2)*X'*y);
V = inv(Sigma0 + (1/Sigma2)*(X'*X));

% posterior draw for coeff|error cov => Gaussian(posterior_mean, posterior_variance)
B = M+(randn(1,1)*chol(V))';
     
end
