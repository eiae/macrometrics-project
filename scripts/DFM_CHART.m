%% DFM charts

% Disclaimer: code written
% by Erik Andres Escayola


%% charts

% plot target along factor in monthly frequency
figure;
subplot(2,1,1);
plot(dates, target, Color=colorDFM, LineWidth=1.2);
axis tight
title('Predicted Real Activity Variable (M)');
subplot(2,1,2);
plot(dates, factor, Color='k', LineWidth=1.2);
hold on
plot(dates,factor_bands(:,[1,end]),':r', LineWidth=0.8);
hold off
axis tight
title('Index of Economic Activity (M)');
sgt = sgtitle('Comparison predicted target variable vs factor of target', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'DFM_predictedTarget_factor.png')); 


% plot target along predicted from common in quarterly frequency
figure;
plot(datesQ, commonQ, Color='r', LineWidth=1.5);
hold on
plot(datesQ, targetQ, Color=colorDFM, LineWidth=1.5)
axis tight
legend('predicted with common (Q)', 'predicted with common & idios (Q)', Location='northwest')
sgt = sgtitle('Predicted target with vs without idiosyncratic component', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'DFM_predictedTarget_AllvsCommon.png')); 


% plot target in mixed frequency
figure;
plot(dates, obs, 'Linewidth',2, 'marker','o', 'color','r');
hold on; 
plot(dates, target, Color=colorDFM, Linewidth=2)
legend('target (Q)', 'predicted target with common & idios (M)', Location='northwest')
axis tight
sgt = sgtitle('Mixed frequency of target vs predicted target', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'DFM_predictedTarget_MixedFreq.png')); 


% plot point and density forecast
figure;
plot(datesForecast, historyPlot, Color=colorDFM, Linewidth=2)
hold on;
plot(datesForecast, forecastPlot, Color=colorDFM, Linewidth=1.5, LineStyle='-')
hold on;
plot(datesForecast(T:end)' , [forecastPlotBands(T:end,1), forecastPlotBands(T:end,end)],...
    Color=colorDFM, LineWidth=1, LineStyle=':')
hFill = [datesForecast(T:end)' fliplr(datesForecast(T:end)')];
inBetween = [forecastPlotBands(T:end,1)', fliplr(forecastPlotBands(T:end,end)')];
fill(hFill , inBetween, colorDFM, FaceAlpha=0.2, LineStyle='none');
legend('history (M)', 'median forecast (M)', sprintf('credible bands %d%%',prct(end)), Location='northwest')
axis tight
grid on
sgt = sgtitle('DFM forecast of predicted target', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'DFM_forecast.png')); 


% check mismatch in scale for observable and predicted target
% figure;
% plot(datesQ, yTargetObs, Color='k', LineWidth=1.5)
% hold on 
% plot(datesQ, yEstimate(nonMissTarget), Color=colorDFM, LineWidth=1.5)
% legend('Quarterly observable target', 'Quarterly predicted target')
% title('Observable vs Predicted')
% 
% disp(corr(yTargetObs,yEstimate(nonMissTarget)))