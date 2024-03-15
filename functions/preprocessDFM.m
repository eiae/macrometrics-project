function [yclean] = preprocessDFM(yraw, selected, T, N, Q)
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
    for i = 6:3:T  % start at Q2 = M6 because of one lag and jump 3 periods 1Q = 3M
        if ~isnan(yraw(i,n))
            % growth rates (period-on-period)
            yclean(i,n) = log(yraw(i,n)/yraw(i-3,n))*100;
        else
            yclean(i,n) = NaN;
        end
    end
    yclean(3,n) = NaN;  % loose first quarterly obs from lag
end

% monthly variables
vars = selected(Q+1:N);
count = 1;
for n = Q+1:N
    check = vars(count);
    if check==6 || check==11 || check==13 || check==14 || check==15 
        yclean(:,n) = yraw(:,n);
    else
        % growth rates (period-on-period)   
        yclean(:,n) = [NaN; diff(log(yraw(:,n)))*100];  % add a NaN since diff() looses one obs
    end
    count = count+1;
end

end
