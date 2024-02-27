function out = lagMaker(x,p)
%% Lag maker
% Input:
% - matrix
% - number of lags
% Output:
% - lagged matrix
% -------------------------------------------------------------------------

% collect size
[R,C] = size(x);

% take first R-p rows of matrix x
x1=x(1:(R-p),:);

% preceed them with p rows of zeros and return
out=[zeros(p,C); x1];

end

