function [B1,B2] = drawCoeffCommon(ztt,Sigma0,B0, lag)
%% Generate draw of autoregressive coeff in common component
% Input:
% - sample of states
% - prior hyperparams 
% Output:
% - draw of autoregressive coeff in common component
% -------------------------------------------------------------------------

% AR estimation for coefficients (Bayesian OLS regression)
Y = ztt;
X = [lagMaker(Y,1) lagMaker(Y,2)];
% remove missing obs
Y = Y(lag+1:end);
X = X(lag+1:end,:);
sigma2 = 1;   % normalized for identification
     
% posterior moments for coeff|erro cov => Gaussian(posterior_mean, posterior_variance)
M = inv(Sigma0+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y);
V = inv(Sigma0+(1/sigma2)*(X'*X));

chck = -1;
while chck<0  % check for stability
    % posterior draw for coeff|error cov => Gaussian(posterior_mean, posterior_variance)
    B = M + (randn(1,2)*chol(V))';
    b = [B(1) B(2); 1 0];
    ee = max(abs(eig(b)));
    if ee<=1
        chck = 1;
    end
end
B1 = B(1,:);
B2 = B(2,:);

end

