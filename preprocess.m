%% charts of original variables

% quarterly charts
varQ = y(:,1:Q);  % get quarterly variables
indexQ = ~isnan(varQ);  % selector for non missing

yQ = varQ(indexQ);  % quarterly variables in quarterly frequency
datesQ = dates(indexQ);  % quarterly dates

nrow = ceil(sqrt(Q));
ncol = ceil(Q / nrow);

figure;
for i = 1:Q
    subplot(nrow, ncol, i)  
    
    plot(datesQ, yQ(:,i),'LineWidth',1.5,'Color','b') 
    
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
    
    plot(dates, y(:,j),'LineWidth',1.5,'Color','k') 
    
    title(names{j}, 'FontSize', 10); 
    axis tight
    grid on
    set(gca,'FontSize',10)
end
sgt = sgtitle('Monthly Series','Interpreter','latex');
sgt.FontSize = 20;


%% transformations


