function [yclean] = preprocess(yraw, dates, names, T, N, Q, M)
%% Preprocess raw data for modelling
% inputs
% - raw dataset
% - dates and names of dataset
% - period length
% - total number of variables
% - number of quarterly and monthly variables
% outputs
% - clean variables after transformations
% -------------------------------------------------------------------------

%% charts of original variables

% quarterly charts
varQ = yraw(:,1:Q);  % get quarterly variables
indexQ = ~isnan(varQ);  % selector for non missing

yrawQ = varQ(indexQ);  % quarterly variables in quarterly frequency
datesQ = dates(indexQ);  % quarterly dates

nrow = ceil(sqrt(Q));
ncol = ceil(Q / nrow);

figure;
for i = 1:Q
    subplot(nrow, ncol, i)  
    
    plot(datesQ, yrawQ(:,i),'LineWidth',1.5,'Color','b') 
    
    title(names{i}, 'FontSize', 10);  
    axis tight
    grid on
    set(gca,'FontSize',10)
end
sgt = sgtitle('Quarterly Series', 'Interpreter','latex');
sgt.FontSize = 20;

% monthly charts
nrow = ceil(sqrt(M));
ncol = ceil(M / nrow);

figure;
for j = Q+1:N  % start iterator at monthly variables
    subplot(nrow, ncol, j-1)  
    
    plot(dates, yraw(:,j),'LineWidth',1.5,'Color','k') 
    
    title(names{j}, 'FontSize', 10); 
    axis tight
    grid on
    set(gca,'FontSize',10)
end
sgt = sgtitle('Monthly Series','Interpreter','latex');
sgt.FontSize = 20;


%% transformations

yclean = yraw;  % copy to have same size
% log levels


% growth rates (period-on-period)
% quarterly variables
for n = 1:Q
    for i = 6:3:T  % start at Q2 = M6 because of one lag and jump 3 periods 1Q = 3M
        if ~isnan(yraw(i,n))
            yclean(i,n) = log(yraw(i,n)/yraw(i-3,n))*100;
        else
            yclean(i,n) = NaN;
        end
    end
    yclean(3,n) = NaN;  % loose first quarterly obs from lag
end



end
















