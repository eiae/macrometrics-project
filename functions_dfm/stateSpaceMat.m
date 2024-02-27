function [Rs,Qs,Hs,Fs] = stateSpaceMat(par, sig2fact, m2q, N)
%% Build state-space param matrices
% Input:
% - model params
% - variance of factor
% - month2quarter conversion vector
% - number of variables
% Output:
% - state-space coeff and innovation/error covariance matrices
% -------------------------------------------------------------------------

% manual construction and preallocation of matrices

% H: coeff matrix obs equation
h01=[m2q'*par(1) m2q' zeros(1,8)];
h2 = ...
[[1 0 0 0 0 0 0 0 ];              
[0 0 1 0 0 0 0 0 ];
[0 0 0 0 1 0 0 0 ];
[0 0 0 0 0 0 1 0 ]];
Hs = [par(2:5) zeros(4,9) h2]; 
Hs = [h01; Hs];  % join cols    

% R: innovation covariance matrix obs equation
Rs = zeros(N,N);                   

% F: coeff matrix state equation
z0 = par(6:7);
z1 = par(8:9);
z2 = par(10:11);
z3 = par(12:13);
z4 = par(14:15);
z5 = par(16:17);  % group model params

f1 = [z0' zeros(1,16)];                
f1a = [eye(4) zeros(4,14)];   
f2 = [zeros(1,5) z1' zeros(1,11)];
f2a = [zeros(4,5) eye(4) zeros(4,9)];
f3 = [zeros(1,10) z2' zeros(1,6)];    
f3a = [zeros(1,10) 1 zeros(1,7)];
f4 = [zeros(1,12) z3' zeros(1,4)];
f4a = [zeros(1,12) 1 zeros(1,5)];
f5 = [zeros(1,14) z4' zeros(1,2)];
f5a = [zeros(1,14) 1 zeros(1,3)];
f6 = [zeros(1,16) z5'];
f6a = [zeros(1,16) 1 zeros(1,1)];  
Fs = [f1;f1a;f2;f2a;f3;f3a;f4;f4a;f5;f5a;f6;f6a];  % join cols

% Q: innovation covariance matrix state equation
q2 = [sig2fact(1) 0 0 0 0 par(18) 0 0 0 0 par(19) 0 par(20) 0 par(21) 0 par(22) 0];     
% q2 = (q2).^2;  % if uncomment run into singularity -> solved with function <invpd()>?
Qs = diag(q2);
   
end
