%% Data preprocessing

% Disclaimer: code written
% by Erik Andres Escayola

% clear
clc


%% original dataset

% quarterly charts
varQ = yraw(:,1:Q);  % get quarterly variables
indexQ = ~isnan(varQ);  % selector for non missing

yrawQ = varQ(indexQ);  % quarterly variables in quarterly frequency
datesQ = dates(indexQ);  % quarterly dates

nrow = ceil(sqrt(Q));
ncol = ceil(Q / nrow);

figure(1);
for i = 1:Q
    subplot(nrow, ncol, i)  
    
    plot(datesQ, yrawQ(:,i),'LineWidth',1.5,'Color','k') 
    
    title(names{i}, 'FontSize', 8);  
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Quarterly Series', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_orig_Q.png')); 

% monthly charts
nrow = ceil(sqrt(M));
ncol = ceil(M / nrow);

figure(2);
for j = Q+1:N  % start iterator at monthly variables
    subplot(nrow, ncol, j-1)  
    
    plot(dates, yraw(:,j),'LineWidth',1.5,'Color','k') 
    
    title(names{j}, 'FontSize', 8); 
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Monthly Series','Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_orig_M.png')); 


%% transformations for DFM

% transformations DFM -> growth rates
yDFM = preprocessDFM(yraw,selected,T,N,Q);

% quarterly charts
varQ = yDFM(:,1:Q);  % get quarterly variables
indexQ = ~isnan(varQ);  % selector for non missing

yDFMQ = varQ(indexQ);  % quarterly variables in quarterly frequency
datesQ = dates(indexQ);  % quarterly dates

nrow = ceil(sqrt(Q));
ncol = ceil(Q / nrow);

figure;
for i = 1:Q
    subplot(nrow, ncol, i)  
    
    plot(datesQ, yDFMQ(:,i),'LineWidth',1.5,'Color',colorDFM) 
    
    title(names{i}, 'FontSize', 8);  
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Quarterly Series (transformed for DFM)', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_DFM_Q.png')); 

% monthly charts
nrow = ceil(sqrt(M));
ncol = ceil(M / nrow);

figure;
for j = Q+1:N  % start iterator at monthly variables
    subplot(nrow, ncol, j-1)  
    
    plot(dates, yDFM(:,j),'LineWidth',1.5,'Color',colorDFM) 
    
    title(names{j}, 'FontSize', 8); 
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Monthly Series (transformed for DFM)','Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_DFM_M.png')); 


%% transformations for VAR

% transformations VAR -> log levels
yVAR = preprocessVAR(yraw,selected,T,N,Q);

% quarterly charts
varQ = yVAR(:,1:Q);  % get quarterly variables
indexQ = ~isnan(varQ);  % selector for non missing

yVARQ = varQ(indexQ);  % quarterly variables in quarterly frequency
datesQ = dates(indexQ);  % quarterly dates

nrow = ceil(sqrt(Q));
ncol = ceil(Q / nrow);

figure;
for i = 1:Q
    subplot(nrow, ncol, i)  
    
    plot(datesQ, yVARQ(:,i),'LineWidth',1.5,'Color',colorVAR) 
    
    title(names{i}, 'FontSize', 8);  
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Quarterly Series (transformed for VAR)', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_VAR_Q.png')); 

% monthly charts
nrow = ceil(sqrt(M));
ncol = ceil(M / nrow);

figure;
for j = Q+1:N  % start iterator at monthly variables
    subplot(nrow, ncol, j-1)  
    
    plot(dates, yVAR(:,j),'LineWidth',1.5,'Color',colorVAR) 
    
    title(names{j}, 'FontSize', 8); 
    axis tight
    grid on
    set(gca,'FontSize',8)
end
sgt = sgtitle('Monthly Series (transformed for VAR)','Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'data_VAR_M.png')); 
