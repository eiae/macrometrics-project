function [sigma2] = drawErrCovIdioLF(ndx,zz,T0,D0,index)
%% Generate draw of error cov of idiosyncratic component for quarterly
% Input:
% - sample of states
% - prior hyperparams 
% - loadings in common component
% - autoregressive coeff idiosyn
% Output:
% - draw of error cov of idiosyncratic component
% -------------------------------------------------------------------------

% use selector vector for idiosyncratic components
z = zz(:,ndx);
Y = delIfCond(z,1-index(:,ndx));

% additional terms for moments
TT = size(Y,1);
      
% posterior moments for error cov|coeff => invGamma(posterior_shape, posterior_scale)
resids = Y;  % residuals
T1 = T0+TT/2;
D1 = (D0+(resids)'*(resids)/2);  % SSR

% posterior draw for error cov|coeff => invGamma(posterior_shape, posterior_scale)
sigma2 = 1/gamrnd(T1,1/D1);  % careful if posterior scale (D1) becomes negative because prior scale (D0) is set to be negative and residuals computation ((resids)'*(resids)/2) too small -> would get NaN draw as scale param should not be negative

end

