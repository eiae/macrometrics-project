function outMat = alternateMat(rows, cols, AR)
%% Build customized matrix with alternating columns
% inputs:
% - number of rows and cols
% - number of autoregressive lags for jumps in cols of state vector (get
% only the contemporaneous elements)
% outputs:
% - customized matrix
% -------------------------------------------------------------------------
outMat = zeros(rows, cols);

for i = 1:rows
    for j = 1:AR:cols
        if  j == i*2-1   % logic of matching 1,2,3,4 rows with 1,3,5,7 cols   
            outMat(i, j) = 1;
        end
    end
end

% example customized matrix 4x8
% [[1 0 0 0 0 0 0 0 ];              
% [0 0 1 0 0 0 0 0 ];
% [0 0 0 0 1 0 0 0 ];
% [0 0 0 0 0 0 1 0 ]];

end