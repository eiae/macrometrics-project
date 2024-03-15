function ans = isposdef(a)
%% Test if matrix positive-definite using Cholesky
% inputs:
% - input matrix
% outputs:
% - ans = 1 is positive definite
% -------------------------------------------------------------------------

[R,p] = chol(a);
ans = (p == 0);

end