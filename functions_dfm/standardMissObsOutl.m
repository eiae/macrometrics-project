function [y,mu,sig] = standardMissObsOutl(x)
%% Standardize that accounting for outliers and missing values
% Input:
% - data matrix with missing values
% Output:
% - standardized vector (accounting for outliers)
% - mean and standard deviation
% -------------------------------------------------------------------------

% specs
[T,n] = size(x);
y = zeros(T,n);  % preallocate

% standardize for each column of matrix (variable)
for i = 1:n
    f = find(1-isnan(x(:,i)));  % filter only non-missing values
    mu(1,i) = mean(x(f,i));
    p25(1,i) = prctile(x(f,i),25);
    p75(1,i) = prctile(x(f,i),75);
    sig(1,i) = std(x(f,i));
    y(f,i) = (x(f,i) - median(x(f,i)))./(p75(1,i) - p25(1,i));  % standardization accounting for outliers 
end

end
