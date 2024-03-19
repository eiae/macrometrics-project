%% Comparison analysis between models

% Disclaimer: code written
% by Erik Andres Escayola

% clear
clc


%% compare RMSE of DFM and VAR

% RMSE values along forecast window
figure;
subplot(2,1,1)
plot(Rmsedates, RmseDFM, Color=colorDFM, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
grid on
axis tight
title('DFM RMSE')

subplot(2,1,2)
plot(RmsedatesAlt, RmseVAR, Color=colorVAR, LineWidth=1.5)
xlabel('forecast window')
ylabel('RMSE')
grid on
axis tight
title('VAR RMSE')
sgt = sgtitle('RMSE values for each forecast iteration', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'RMSE_values.png')); 

% RMSE kernel densities
figure;
subplot(2,1,1)
plot(xiRmseDFM, fRmseDFM, Color=colorDFM, LineWidth=1.5)
xlabel('RMSE')
ylabel('density')
title('DFM RMSE kernel')

subplot(2,1,2)
plot(xiRmseVAR, fRmseVAR, Color=colorVAR, LineWidth=1.5)
xlabel('RMSE')
ylabel('density')
title('VAR RMSE kernel')
sgt = sgtitle('RMSE kernels for total forecast iterations', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'RMSE_kernel.png')); 
