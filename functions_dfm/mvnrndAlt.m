function r = mvnrndAlt(mu,sigma,cases)
%% Modified matlab function to be able to draw from almost singular matrices
% inputs
% - mean, variance
% - number of rows in matrix
% outputs
% - draws from multivariate random normal
% -------------------------------------------------------------------------

% check if objects ill-defined
[m1 n1] = size(mu);
c = max([m1 n1]);
if m1 .* n1 ~= c
   error('Mu must be a vector.');
end

[m n] = size(sigma);
if m ~= n
   error('Sigma must be square');
end

if m ~= c
   error('The length of mu must equal the number of rows in sigma.');
end


[T p] = chol(sigma);
if p > 0
   % The input covariance has some type of perfect correlation.
   % If it is positive semi-definite, we can still proceed.
   % Find a factor T that has the property    sigma == T' * T.
   % Unlike a Cholesky factor, this T is not necessarily triangular
   % or even square.
   T = mvnfactor(sigma);
   %if p > 0, error('Sigma must be a non-negative definite matrix.'); end
end

if m1 == c
  mu = mu';
end

mu = mu(ones(cases,1),:);

r = randn(cases,size(T,1)) * T + mu;

% -------------------------------------------------------------------------

function T = mvnfactor(sigma) % originally [T,p] = 
%MVNFACTOR  Do Cholesky-like decomposition, allowing zero eigenvalues
%   SIGMA must be symmetric.  In general T is not square or triangular.
[U,D] = eig((sigma+sigma')/2);

D = diag(D);
tol = max(D) * length(D) * eps;
t = (D > tol);
D = D(t);
T = diag(sqrt(D)) * U(:,t)';

end
end