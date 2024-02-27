function y = delIfCond(x,cond)
%% Select values for which condition is false
% Input:
% - input vec
% - condition as bolean vector 
% Output:
% - output vec
% -------------------------------------------------------------------------
       
% select values if condition false
y = x(cond==0);

end