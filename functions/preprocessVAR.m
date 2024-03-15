function [yclean] = preprocessVAR(yraw, selected, T, N, Q)
%% Preprocess raw data for modelling
% inputs
% - raw dataset
% - period length
% - total number of variables
% - number of quarterly variables
% outputs
% - clean variables after transformations
% -------------------------------------------------------------------------

yclean = yraw;  % copy to have same size

% quarterly variables
for n = 1:Q
    for i = 3:3:T  % start at Q1 = M3 because of jump 3 periods 1Q = 3M
        if ~isnan(yraw(i,n))
            % log
            yclean(i,n) = log(yraw(i,n));
        else
            yclean(i,n) = NaN;
        end
    end
end

% monthly variables
vars = selected(Q+1:N);
count = 1;
for n = Q+1:N
    check = vars(count);
    if check==6 || check==9 || check==13 
        yclean(:,n) = yraw(:,n);
    else
        % log 
        yclean(:,n) = log(yraw(:,n));  
    end
    count = count+1;
end

end
