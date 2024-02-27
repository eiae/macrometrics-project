function  [xinv,flag] = invNonPdMat(x)
%% Inverts matrix using augmented eigenvalues for non positive-definite case
% inputs:
% - input matrix
% outputs:
% - flag (=0 pd, =1 non pd, eigenval augmented
% - inverted matrix
% -------------------------------------------------------------------------

if isposdef(x)
    xinv = inv(x);
    flag = 0;
else
    % carries out SVD and augments small eigenvalues, setting them to 0.01
    % this is intended to ensure PD var-cov matrices returned by numerical hessians
    flag = 1;
    [n m] = size(x);
    [u d v] = svd(x);
    dd = diag(d);
    xchk = u*diag(dd)*v';
    dd = dd + 1000*eps;
    di = ones(n,1)./dd;
    xinv = u*diag(di)*v';
end

end

