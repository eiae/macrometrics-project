%% VAR charts

% Disclaimer: code written
% by Erik Andres Escayola


%% charts 

% plot target in mixed frequency
figure;
plot(dates, obsAlt, 'Linewidth',2, 'marker','o', 'color','r');
hold on; 
plot(dates, targetAlt, 'Linewidth',2, 'color', colorVAR)
legend('target (Q)', 'predicted target (M)', Location='northwest')
axis tight
sgt = sgtitle('Mixed frequency of target vs predicted target', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'VAR_predictedTarget_MixedFreq.png')); 


% plot point and density forecast
figure;
plot(datesForecast, historyPlotAlt, Color=colorVAR, Linewidth=2)
hold on;
plot(datesForecast, forecastPlotAlt, Color=colorVAR, Linewidth=1.5, LineStyle='-')
hold on;
plot(datesForecast(T:end)' , [forecastPlotBandsAlt(T:end,1), forecastPlotBandsAlt(T:end,end)],...
    Color=colorVAR, LineWidth=1, LineStyle=':')
hFill = [datesForecast(T:end)' fliplr(datesForecast(T:end)')];
inBetween = [forecastPlotBandsAlt(T:end,1)', fliplr(forecastPlotBandsAlt(T:end,end)')];
fill(hFill , inBetween, colorVAR, FaceAlpha=0.2, LineStyle='none');
legend('history (M)', 'median forecast (M)', sprintf('credible bands %d%%',prct(end)), Location='northwest')
axis tight
grid on
sgt = sgtitle('VAR forecast of predicted target', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'VAR_forecast.png')); 


% levels
% figure;
% plot(dates, obsAltLvl, 'Linewidth',2, 'marker','o', 'color','r'); 
% hold on; 
% plot(dates, targetAltLvl, 'Linewidth',2, 'color', colorVAR)
% legend('observed (quarterly)', 'index (monthly)')
% axis tight
% title('Quarterly Level of Observed vs Index of Economic Activity');
