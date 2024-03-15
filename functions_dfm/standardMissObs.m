function [y,mu,sig] = standardMissObs(x)
%% Standardize that accounting for missing values
% Input:
% - data matrix with missing values
% Output:
% - standardized vector 
% - mean and standard deviation
% -------------------------------------------------------------------------

% specs
[T,n] = size(x);
y = zeros(T,n);  % preallocate

% standardize for each column of matrix (variable)
for i = 1:n
    f = find(1-isnan(x(:,i)));  % filter only non-missing values
    mu(1,i) = mean(x(f,i));
    sig(1,i) = std(x(f,i));
    y(f,i) = (x(f,i) - mu(1,i))/sig(1,i);  % standardization 
end

end