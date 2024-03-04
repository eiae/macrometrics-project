%% MF Bayesian VAR with unbalanced data

% Disclaimer: code adapted from Empirical Macro Toolbox
% by Filippo Ferroni and Fabio Canova
% https://github.com/naffe15/BVAR_/tree/master

close all
clc
seed=0;  
rng(seed);  

 
% load the mixed frequency data
%load DataMF
%y = [GDP IPI HICP CORE Euribor1Y UNRATE]; % variables  to  be  used
lags = 2;                                 % number  of  lags 

% T aggregation: the quarterly variable  
% xq(t) = 1/3( xm(t) + xm(t-1) + xm(t-2)) at least two lags are needed
options.mf_varindex     = 1; 
options.K               = 1000;   
options.priors.name     = 'Minnesota';
options.noprint         = 1;
bvarmf                  = bvar_(yVAR,lags,options); % estimate the var model

sorty = sort(bvarmf.yfill,3);   % sort the smoothed  state


% plot: levels
figure('Name','Monthly EA GDP')
plot(dates,yVAR(:,1),'Linewidth',2,'marker','o','color','r'); hold on; 
plot(dates,sorty(:,1,0.5*bvarmf.ndraws),'b','Linewidth',2)
legend('Quarterly data','MXFREQ VAR flow')
axis tight

% plot: growth rates
Dy = 100*(yVAR(4:end,1)-yVAR(1:end-3,1));
Y = sorty(:,1,0.5*bvarmf.ndraws);
DY = 100*(Y(4:end,1)-Y(1:end-3,1));

figure('Name','Monthly EA GDP growth')
plot(dates(4:end),Dy,'ro');
hold on; plot(dates(4:end),DY,'b','Linewidth',2)
legend('Quarterly data','MXFREQ VAR flow ')
axis tight
