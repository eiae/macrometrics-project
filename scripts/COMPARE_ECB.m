%% Comparison analysis with ECB projections

% Disclaimer: code written
% by Erik Andres Escayola

% clear
clc


%% ECB projections
figure;
plot(datesECB(end-periodsBack:end), dataECB(end-periodsBack:end,2), Color=[1 0.71 0], LineWidth=1.5)
hold on
plot(datesECB(end-periodsBack:end), dataECB(end-periodsBack:end,1), Color=[0 0.22 0.6], LineWidth=1.5)
axis tight
grid on 
legend(namesECB{2}, namesECB{1}, Location='northeast')
sgt = sgtitle('ECB projections on euro area real GDP growth', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'ECB_projections.png')); 


%% compare DFM and ECB
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
plot(datesECBCompare(end-periodsBack:end), dataECBCompare(end-periodsBack:end,2), Color=[1 0.71 0], LineWidth=1.5)
hold on
plot(datesECBCompare(end-periodsBack:end), dataECBCompare(end-periodsBack:end,1), Color=[0 0.22 0.6], LineWidth=1.5)
plot(datesECBCompare(end-periodsBack:end), forecastQ(end-periodsBack:end), Color=colorDFM, LineWidth=1.5)
hFill = [datesECBCompare(end-periodsBack:end)' fliplr(datesECBCompare(end-periodsBack:end)')];
inBetween = [forecast_bandsQ(end-periodsBack:end,1)', fliplr(forecast_bandsQ(end-periodsBack:end,end)')];
fill(hFill , inBetween, colorDFM, FaceAlpha=0.2, LineStyle='none');
line([datesQ(end), datesQ(end)], ylim, Color='k', LineStyle='--', LineWidth=1);  % add vertical line when forecast starts (last quarter of data)
axis tight
grid on 
legend(namesECB{2}, namesECB{1}, 'median forecast DFM', sprintf('credible bands %d%%',prct(end)), Location='northwest')
sgt = sgtitle('Comparison ECB projections and DFM forecast', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'Compare_ECB_DFM_forecast.png')); 


%% compare VAR to ECB
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
plot(datesECBCompare(end-periodsBack:end), dataECBCompare(end-periodsBack:end,2), Color=[1 0.71 0], LineWidth=1.5)
hold on
plot(datesECBCompare(end-periodsBack:end), dataECBCompare(end-periodsBack:end,1), Color=[0 0.22 0.6], LineWidth=1.5)
plot(datesECBCompare(end-periodsBack:end), forecastAltQ(end-periodsBack:end), Color=colorVAR, LineWidth=1.5)
hFill = [datesECBCompare(end-periodsBack:end)' fliplr(datesECBCompare(end-periodsBack:end)')];
inBetween = [forecast_bandsAltQ(end-periodsBack:end,1)', fliplr(forecast_bandsAltQ(end-periodsBack:end,end)')];
fill(hFill , inBetween, colorVAR, FaceAlpha=0.2, LineStyle='none');
line([datesQ(end), datesQ(end)], ylim, Color='k', LineStyle='--', LineWidth=1);  % add vertical line when forecast starts (last quarter of data)
axis tight
grid on 
legend(namesECB{2}, namesECB{1}, 'median forecast VAR', sprintf('credible bands %d%%',prct(end)), Location='northwest')
sgt = sgtitle('Comparison ECB projections and VAR forecast', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'Compare_ECB_VAR_forecast.png')); 


%% compare DFM and VAR to ECB

% compute deltas and correlations
forecastStart = length(datesECB) - length(datesECBCompare(end-horizonECB:end));

% delta -> difference between ECB (Mar24-MPE) and models forecast for horizons
forecastDeltaECBvsDFM = dataECBCompare(forecastStart:end,1) - forecastQ(forecastStart:end);  % DFM
[fECBvsDFM, xiECBvsDFM] = ksdensity(forecastDeltaECBvsDFM);
forecastDeltaECBvsVAR = dataECBCompare(forecastStart:end,1) - forecastAltQ(forecastStart:end);  % VAR
[fECBvsVAR, xiECBvsVAR] = ksdensity(forecastDeltaECBvsVAR);

figure;
plot(xiECBvsDFM,fECBvsDFM, Color=colorDFM, LineWidth=1.5);
hold on
plot(xiECBvsVAR,fECBvsVAR, Color=colorVAR, LineWidth=1.5);
legend('DFM forecast difference density', 'VAR forecast difference density ', Location='northwest')
sgt = sgtitle('Kernel density of forecast differences: ECB Mar23 MPE vs models', 'Interpreter','latex');
sgt.FontSize = 12;
saveas(gcf, fullfile(savepath,'Compare_ECB_DFM_VAR_errorKernel.png')); 

% correlation -> relation between ECB (Mar24-MPE) and DFM forecast for horizons
forecastCorrECBvsDFM = corr(dataECBCompare(forecastStart:end,1), forecastQ(forecastStart:end));
fprintf('Forecast correlation between ECB Mar23 MPE and DFM:\n %4.2f \n', forecastCorrECBvsDFM)

% correlation -> relation between ECB (Mar24-MPE) and VAR forecast for horizons
forecastCorrECBvsVAR = corr(dataECBCompare(forecastStart:end,1), forecastAltQ(forecastStart:end));
fprintf('Forecast correlation between ECB Mar23 MPE and VAR: \n %4.2f \n', forecastCorrECBvsVAR)