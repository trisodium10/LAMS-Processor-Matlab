
% Plot to compare LAMS (points), Kalman (line), and raw radome (triangles)
figure; 
plot(vLAMSac.','.');
set(gca, 'ColorOrderIndex',1);
hold on;
plot(vRD.','^');
set(gca, 'ColorOrderIndex',1);
plot(vRDEst.','-','linewidth',1.5);
title(['Velocity Components LAMS vs. RD, ', ncfilenameB], 'Interpreter','none');
grid on;
xlabel('Index');
ylabel('Velocity [m/s]');
legend('x','y','z')

% Plot to compare LAMS (blue points), Kalman (yellow dashed), and raw radome (red solid)
% with each component displayed on separate plots
figure; 
subplot(3,1,1)
plot(timeWL/3600,vLAMSac(1,:).*prod(foundNan,1),'.');
hold on;
plot(timeWL/3600,vRD(1,:),'-','linewidth',1.5);
plot(timeWL/3600,vRDEst(1,:),'--','linewidth',1.5);
grid on;
ylabel('v_x [m/s]');
legend('LAMS','Radome','Kalman');

subplot(3,1,2)
plot(timeWL/3600,vLAMSac(2,:).*prod(foundNan,1),'.');
hold on;
plot(timeWL/3600,vRD(2,:),'-','linewidth',1.5);
plot(timeWL/3600,vRDEst(2,:),'--','linewidth',1.5);
grid on;
ylabel('v_y [m/s]');

subplot(3,1,3)
plot(timeWL/3600,vLAMSac(3,:).*prod(foundNan,1),'.');
hold on;
plot(timeWL/3600,vRD(3,:),'-','linewidth',1.5);
plot(timeWL/3600,vRDEst(3,:),'--','linewidth',1.5);
grid on;
ylabel('v_z [m/s]');
xlabel('Time [h-UTC]');

set(gcf,'NextPlot','add');
axes;
h = title(ncfilenameB(1:end-3),'Interpreter','none');
set(gca,'Visible','off');
set(h,'Visible','on');



% figure; 
% plot(EstErrors.')
% title('LAMS-RD Error');
% title(['LAMS-RD Error, ', ncfilenameB], 'Interpreter','none');
% 
% figure; 
% plot((aDiff.*foundPk).')
% title('LAMS-RD LOS Error');
% title(['LAMS-RD LOS Error, ', ncfilenameB], 'Interpreter','none');

% Displays difference between LAMS-Kalman Filtered Radome Beam LOS spees 
% where LAMS data exists.
figure; 
plot(timeWL/3600,(aDiff.*foundNan).','linewidth',1.5)
xlabel('Time [h-UTC]');
ylabel('Speed Difference [m/s]');
grid on;
legend('Beam 1','Beam 2','Beam 3');
title(['LAMS-RD LOS Speed Difference, ', ncfilenameB], 'Interpreter','none');

% Plot time series of sensitivity coefficients
figure;
plot(timeWL/3600,CoeffRecord.');
xlabel('Time [h-UTC]');
ylabel('Coefficient Value');
grid on;
legend(CoefLabel);
title(strcat("Aircraft Coefficients, ", ncfilenameB),'Interpreter','none');

figure;
for ai = 1:size(CoeffRecord,1)
    bincoef = linspace(min(CoeffRecord(ai,:)),max(CoeffRecord(ai,:)),40);
    dbin = mean(diff(bincoef));
    bincoef = [bincoef, dbin+bincoef(end)] - dbin/2.0;
    hcoef = hist(CoeffRecord(ai,:),bincoef);
    plot(bincoef,hcoef);
    hold on;
end
grid on;
title(strcat("Aircraft Coefficients, ", ncfilenameB),'Interpreter','none');
legend(CoefLabel);

% Shows the estimated uncertainty of the x,y,z airflow components in the 
% LAMS coordinate frame.  Includes Kalman estimated error components.
figure;
semilogy(timeWL/3600,sqrt(covLAMSwl.'));
legend('\sigma_{vx}','\sigma_{vy}','\sigma_{vz}')
xlabel('Time [h-UTC]');
ylabel('Velocity Uncertainty [m/s]')
grid on;
title(['Estimated 3 Component Uncertainty in ' ncfilenameB],'Interpreter','none');

% Shows the estimated uncertainty of the beam LOS speeds
% Includes Kalman estimated error components.
figure;
semilogy(timeWL/3600,sqrt(covLAMSLOS.'));
legend('\sigma_{vx}','\sigma_{vy}','\sigma_{vz}')
xlabel('Time [h-UTC]');
ylabel('Velocity Uncertainty [m/s]')
grid on;
title(['Estimated LOS Uncertainty in ' ncfilenameB],'Interpreter','none');


% Plots the statistics on the difference between LOS LAMS observations
% and Kalman predicted LOS peak locations.
% Also reports the fraction of data points where a LAMS peak was found
PkFoundStats = sum(foundPk,2)./size(foundPk,2);
hbins = -3:0.01:3;
hDifflos = zeros(length(hbins),sum(BeamList));
PltStr = '';
for ai = 1:sum(BeamList)
    DiffData = aDiff(ai,:);
    DiffData(foundPk(ai,:)==0) = nan;
    hDifflos(:,ai) = hist(DiffData,hbins);
    legStr{ai} = strcat("Beam ", num2str(ai));
    PltStr = strcat(PltStr, "\sigma_", num2str(ai), " = ", num2str(nanstd(DiffData)), " m/s \newline");
    PltStr = strcat(PltStr, "Beam ", num2str(ai), " Peaks Found: ", num2str(PkFoundStats(ai)*100), "%\newline");
end
figure; 
plot(hbins,hDifflos.','linewidth',2);
xlabel('LOS Kalman Error [m/s]');
ylabel('Count')
grid on;
legend(legStr);
text(hbins(1),ylim()*[0;0.8],PltStr);
title(ncfilenameB,'Interpreter','none');

% Compare pitot-static TASX
%         TAS from the magnitude of the LAMS velocity vector
%         TAS from the magnitude of the Kalman velocity vector
% Also plots error envelope around LAMS TAS
figure; 
plot(timeWL/3600,[tas(iD0-1+(1:DataLen)) tasLAMS.' tasEstv.'])
hold on;
plot(timeWL/3600,(tasLAMS.')+[-sqrt(var_tasLAMS.'),sqrt(var_tasLAMS.')],'k--'),
xlabel('Time [h UT]');
ylabel('TAS, |v| [m/s]');
legend('TASX','LAMS (magnitude)','Kalman (magnitude)','std LAMS/Kalman');
grid on;
title(ncfilenameB,'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_TAS_Time', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_TAS_Time', FileTag, '.png'));
end


% histogram of difference between TAS and LAMS
Diff_Edges = -4:0.01:4;

herr2 = hist(errorTotalAirspeed2.*prod(foundNan,1),Diff_Edges);

figure;
plot(Diff_Edges,herr2/sum(herr2),'linewidth',2);
ylabel('Frequency');
xlabel('Difference in True Airspeed [m/s]');
title(ncfilenameB,'Interpreter','none');
grid on;
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TAS', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TAS', FileTag, '.png'));
end


% Plots the difference in TAS between LAMS and TASX and LAMS and Kalman
figure; 
plot(timeWL/3600,[sqrt(sum((vRDEst).^2,1)).'-sqrt(sum((vLAMSwl).^2,1)).',(errorTotalAirspeed.*prod(foundNan,1)).'])
xlabel('Mission Time [h UT]');
ylabel('Difference in Total Airspeed [m/s]');
legend('Kalman-LAMS','RD-LAMS','TAS-LAMS'); %,'std LAMS/Kalman');
title(ncfilenameB,'Interpreter','none');
grid on;
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TotalAirspeed_Time', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TotalAirspeed_Time', FileTag, '.png'));
end

% % Comparisons for different methods computing TAS.  Generally TAS is
% % assumed to be the magnitude of the full vector of airflow and it is
% % assumed that the pitot-static system provides this without factoring in
% % sideslip and angle of attack components.  However I found I got better
% % agreement with LAMS if I factored in those components separately and
% % treated TASX as only the foward velocity component.  At this stage, the
% % processing treats TASX as the full vector and the radome observed vector
% % is normalized to equal TASX.
% figure; 
% plot(timeWL/3600,[(errorTotalAirspeed.*prod(foundNan,1)).',(errorTotalAirspeed2.*prod(foundNan,1)).',(errorTotalAirspeed3.*foundNan(2,:)).',(errorTotalAirspeed4.*foundNan(2,:)).'])
% xlabel('Mission Time [h UT]');
% ylabel('Difference in True Airspeed [m/s]');
% legend('RD-LAMS','TAS-LAMS','RD-LAMS Forward','TAS-LAMS Forward'); %,'std LAMS/Kalman');
% title(ncfilenameB,'Interpreter','none');
% grid on;
% 
% Diff_Edges = -4:0.01:4;
% herr1 = hist(errorTotalAirspeed.*prod(foundNan,1),Diff_Edges);
% herr2 = hist(errorTotalAirspeed2.*prod(foundNan,1),Diff_Edges);
% herr3 = hist(errorTotalAirspeed3.*foundNan(2,:),Diff_Edges);
% herr4 = hist(errorTotalAirspeed4.*foundNan(2,:),Diff_Edges);
% 
% figure;
% plot(Diff_Edges,herr1/sum(herr1)); hold on;
% plot(Diff_Edges,herr2/sum(herr2));
% plot(Diff_Edges,herr3/sum(herr3));
% plot(Diff_Edges,herr4/sum(herr4));
% ylabel('Frequency');
% xlabel('Difference in Total Airspeed [m/s]');
% legend('RD-LAMS','TAS-LAMS','RD-LAMS Forward','TAS-LAMS Forward'); %,'std LAMS/Kalman');
% title(ncfilenameB,'Interpreter','none');
% grid on;

% Combined plot of TASX, Kalman TAS and LAMS TAS (top plot)
% with differences with LAMS (bottom plot)
figure; 
subplot(2,1,1);
plot(timeWL/3600,sqrt(sum(vRDEst.^2,1)).','-'); hold on;
plot(timeWL/3600,sqrt(sum(vRD.^2,1)).');
plot(timeWL/3600,sqrt(sum(vLAMSwl.^2,1)).','--');
ylabel('|v| [m/s]');
legend('Kalman','Radome','LAMS'); %,'std LAMS/Kalman');
grid on;
subplot(2,1,2);
plot(timeWL/3600,[sqrt(sum(vRDEst.^2,1).*prod(foundNan,1)).'-sqrt(sum(vLAMSwl.^2,1)).',(errorTotalAirspeed.*prod(foundNan,1)).'])
xlabel('Mission Time [h UT]');
ylabel('\Delta |v| [m/s]');
legend('Kalman-LAMS','Radome-LAMS');
grid on;
ylim([-5,5]);

% Compare the angle difference between Radome and LAMS and Kalman and LAMS
% Notably, the radome tends to do better than Kalman
windAngleDiff = acos(sum((vLAMSac./sqrt(sum(vLAMSac.^2,1))).*(vRD./sqrt(sum(vRD.^2,1))),1)).*prod(foundNan,1);
windAngleDiffEst = acos(sum((vLAMSac./sqrt(sum(vLAMSac.^2,1))).*(vRD./sqrt(sum(vRDEst.^2,1))),1)).*prod(foundNan,1);
figure; 
plot(timeWL/3600,[windAngleDiff.' windAngleDiffEst.']*180/pi);
grid on;
xlabel('Mission Time [h UT]');
ylabel('Angle Between Vectors [deg.]');
legend('RD-LAMS','Kalman-LAMS');
title(ncfilenameB,'Interpreter','none');

% Scatter plot to look for correlations in errors of TAS and differences
% in vector direction between the observations
figure; 
scatter(windAngleDiff*180/pi,errorTotalAirspeed); 
hold on; 
scatter(windAngleDiffEst*180/pi,sqrt(sum(vRDEst.^2,1))-sqrt(sum(vLAMSwl.^2,1)));
grid on;
xlabel('Angle Between Vectors [deg.]');
ylabel('Difference in Total Airspeed [m/s]');
title(ncfilenameB,'Interpreter','none');



EdgeWinds = -15:0.05:15;
% Create histograms to compare measurements of vertical winds
hwic1 = hist(wic(iD0-1+(1:DataLen)).',EdgeWinds); 
hwicL = hist(wLAMS(3,:).*prod(foundNan,1),EdgeWinds);
hwicK = hist(wEst(3,:),EdgeWinds);
figure; 
plot(EdgeWinds,[hwic1; hwicL; hwicK].','linewidth',2);
grid on;
grid minor
legend('WIC','LAMS','Kalman');
xlabel('Vertical Wind [m/s]'); 
ylabel('Counts');
xlim([-3,3])
title(['Vertical Winds, ' ncfilenameB],'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Vert', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Vert', FileTag, '.png'));
end

% % Time series of vertical winds
% figure;
% plot(timeWL/3600,[wic(iD0-1+(1:DataLen)).'; wLAMS(3,:); wEst(3,:)].','linewidth',1.5);
% hold on;
% plot(timeWL/3600,0.5*sum([wLAMS(3,:); wEst(3,:)].',2)+sqrt(covLAMSglob(3,:).').*[-1, 1],'k--');
% legend('WIC','LAMS','Kalman','std LAMS/Kalman');
% ylim([-10 10]);
% xlabel('Mission Time [h UT]');
% title(['Vertical Winds, ' ncfilenameB],'Interpreter','none');
% if SaveFigs
%    saveas(gcf,[ncfilenameB(1:end-3) '_LAMS_Wind_Vert_Time' FileTag '.fig']);
%    saveas(gcf,[ncfilenameB(1:end-3) '_LAMS_Wind_Vert_Time' FileTag '.png']);
% end


% Aircraft frame comparison of three component winds between Kalman and
% LAMS
hbins = -3:0.01:3;
hDiff = zeros(length(hbins),3);
CoordList = {'x','y','z'};
PltStr = '';
Diffstd = zeros(3,1);
for ai = 1:3
    DiffData = vLAMSac(ai,:) - vRDEst(ai,:);
    DiffData(prod(foundPk,1)==0) = nan;
    hDiff(:,ai) = hist(DiffData,hbins);
    legStr{ai} = [CoordList{ai}];
    Diffstd(ai) = nanstd(DiffData);
    PltStr = strcat(PltStr, '\sigma_', CoordList{ai}, ' = ', num2str(nanstd(DiffData)), ' m/s \newline');
end
figure; 
plot(hbins,hDiff.','linewidth',2);
xlabel('Three Component Kalman Error [m/s]');
ylabel('Count')
grid on;
legend(legStr);
text(hbins(1)*0.9,ylim()*[0;0.8],PltStr);
title(ncfilenameB,'Interpreter','none');

% Aircraft frame comparison of three component winds between Kalman and
% LAMS AND Radome and LAMS
PkFoundStats = sum(foundPk,2)./size(foundPk,2);
hbins = -3:0.01:3;
hDiff2 = zeros(length(hbins),sum(BeamList));
for ai = 1:sum(BeamList)
    DiffData = vLAMSac(ai,:) - vRD(ai,:);
    DiffData(prod(foundPk,1)==0) = nan;
    hDiff2(:,ai) = hist(DiffData,hbins);
    legStr{ai+3} = ['Radome ' CoordList{ai}];
    PltStr = strcat(PltStr, 'Radome \sigma_', CoordList{ai}, ' = ', num2str(nanstd(DiffData)), ' m/s \newline');
end
figure; 
plot(hbins,hDiff.','linewidth',2);
set(gca, 'ColorOrderIndex',1);
hold on;
plot(hbins,hDiff2.','--','linewidth',2);
xlabel('Three Component Airspeed Difference [m/s]');
ylabel('Count')
grid on;
grid minor;
legend(legStr);
text(hbins(1),ylim()*[0;0.8],PltStr);
title(ncfilenameB,'Interpreter','none');



% radome estimated wind vector 
wRD = [(wsc(iD0-1+(1:DataLen)).*cos(wdc(iD0-1+(1:DataLen)))).'; ...
    (wsc(iD0-1+(1:DataLen)).*sin(wdc(iD0-1+(1:DataLen)))).'; ...
    wic(iD0-1+(1:DataLen)).'];

EdgeWinds = -15:0.05:15;
% Create histograms to compare measurements of vertical winds
hwic1 = hist((wRD(3,:)-wLAMS(3,:).*prod(foundNan,1)).',EdgeWinds); 
hwNS = hist((wRD(1,:)-wLAMS(1,:).*prod(foundNan,1)).',EdgeWinds);
hwEW = hist((wRD(2,:)-wLAMS(2,:).*prod(foundNan,1)).',EdgeWinds);
figure; 
plot(EdgeWinds,[hwNS; hwEW; hwic1].','linewidth',2);
grid on;
grid minor
legend('N/S','E/W','z');
xlabel('Difference in Winds [m/s]'); 
ylabel('Counts');
xlim([-3,3])
title(['Radome - LAMS Winds, ' ncfilenameB],'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff', FileTag, '.png'));
end



% radome estimated wind vector of sub domain
% [~,i1] = min(abs(timeWL-3600*22.1)); % rf01
% [~,i2] = min(abs(timeWL-3600*23.1));
% [~,i1] = min(abs(timeWL-3600*21.6)); % rf03
% [~,i2] = min(abs(timeWL-3600*23));
% [~,i1] = min(abs(timeWL-3600*0));
% [~,i2] = min(abs(timeWL-3600*50));

% [~,i1] = min(abs(timeWL-3600*22-7*60)); % rf03
% [~,i2] = min(abs(timeWL-3600*22-7*60-180));
% [~,i1] = min(abs(timeWL-3600*22-16*60-40)); % rf03
% [~,i2] = min(abs(timeWL-3600*22-16*60-40-180));
% [~,i1] = min(abs(timeWL-3600*22-3*60)); % rf03
% [~,i2] = min(abs(timeWL-3600*22-3*60-180));
[~,i1] = min(abs(timeWL-3600*22-10*60-30)); % rf03
[~,i2] = min(abs(timeWL-3600*22-10*60-30-180));


wdc2 = wdc(iD0-1+(1:DataLen));
thdg2 = thdg(iD0-1+(1:DataLen));

EdgeWinds = -15:0.05:15;
% Create histograms to compare measurements of vertical winds
hwic1 = hist((wRD(3,i1:i2)-wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds); 
hwNS = hist((wRD(1,i1:i2)-wLAMS(1,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds);
hwEW = hist((wRD(2,i1:i2)-wLAMS(2,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds);
figure; 
plot(EdgeWinds,[hwNS; hwEW; hwic1].','linewidth',2);
grid on;
grid minor
legend('N/S','E/W','z');
xlabel('Difference in Winds [m/s]'); 
ylabel('Counts');
xlim([-3,3])
title(['Radome - LAMS Winds, ' ncfilenameB],'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.png'));
end


EdgeWinds = -15:0.05:15;
% Create histograms to compare measurements of vertical winds
hwic1 = hist((wRD(3,i1:i2)-wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds); 
hwHoriz = hist(((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)).',EdgeWinds);

figure; 
plot(EdgeWinds,[hwHoriz;hwic1].','linewidth',2);
grid on;
grid minor
legend('horizontal','vertical');
xlabel('Difference in Winds [m/s]'); 
ylabel('Counts');
xlim([-3,3])
title(['Radome - LAMS Winds, ' ncfilenameB],'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff_hv_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Wind_Diff_hv_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.png'));
end

hWind2d = histogram2((wRD(1,i1:i2)-wLAMS(1,i1:i2).*prod(foundNan(:,i1:i2),1)).',(wRD(2,i1:i2)-wLAMS(2,i1:i2).*prod(foundNan(:,i1:i2),1)).','XBinEdges',EdgeWinds,'YBinEdges',EdgeWinds);


hwicLAMS = hist((wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds); 
hwicRD = hist((wRD(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',EdgeWinds); 
figure; 
plot(EdgeWinds,[hwicRD ;hwicLAMS].','linewidth',2);
grid on;
grid minor
legend('Radome','LAMS');
xlabel('Vertical Wind [m/s]'); 
ylabel('Counts');
xlim([-3,3])
title(['Vertical Winds, ' ncfilenameB],'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_VerticalWind_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_VerticalWind_', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.png'));
end

% histogram of difference between TAS and LAMS
Diff_Edges = -4:0.01:4;

herr2 = hist(errorTotalAirspeed2(i1:i2).*prod(foundNan(:,i1:i2),1),Diff_Edges);

figure;
plot(Diff_Edges,herr2/sum(herr2),'linewidth',2);
ylabel('Frequency');
xlabel('Difference in True Airspeed [m/s]');
title(ncfilenameB,'Interpreter','none');
grid on;
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TAS', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Diff_TAS', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.png'));
end


% figure; 
% scatter(sqrt(sum(wRD(1:2,i1:i2).^2,1).').*(cos(wdc2(i1:i2))-cos(thdg2(i1:i2))),sqrt(sum(wRD(1:2,i1:i2).^2,1).').*(sin(wdc2(i1:i2))-sin(thdg2(i1:i2))),1,((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)))

% scatter(cos(wdc2(i1:i2))-cos(thdg2(i1:i2)),sin(wdc2(i1:i2))-sin(thdg2(i1:i2)),sqrt(sum(wRD(1:2,i1:i2).^2,1)),((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)))
% scatter(cos(wdc2(i1:i2)-thdg2(i1:i2)),sin(wdc2(i1:i2)-thdg2(i1:i2)),sqrt(sum(wRD(1:2,i1:i2).^2,1)),((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)))
figure; 
scatter(sqrt(sum(wRD(1:2,i1:i2).^2,1).').*cos(wdc2(i1:i2)-thdg2(i1:i2)),sqrt(sum(wRD(1:2,i1:i2).^2,1).').*sin(wdc2(i1:i2)-thdg2(i1:i2)),sqrt(sum(wRD(1:2,i1:i2).^2,1)),((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)))
caxis([-1,1])
grid on;
xlabel('Tailwind [m/s]');
ylabel('Cross wind [m/s]');
colorbar;
axis equal
title(ncfilenameB,'Interpreter','none');

figure; 
scatter(wRD(1,i1:i2),wRD(2,i1:i2),sqrt(sum(wRD(1:2,i1:i2).^2,1)),((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)))
caxis([-1,1])
grid on;
xlabel('Tailwind [m/s]');
ylabel('Cross wind [m/s]');
colorbar;
axis equal
title(ncfilenameB,'Interpreter','none');

figure; 
scatter(vRD(1,:),vRD(2,:),1,(tasLAMS.'-tas(iD0-1+(1:DataLen))).*prod(foundNan,1).')
caxis([-1,1])
grid on;
xlabel('Tail airflow [m/s]');
ylabel('Cross airflow [m/s]');
colorbar;
% axis equal
title(ncfilenameB,'Interpreter','none');

hWbins = linspace(-3,3,200);
hwerr = hist3([((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)).',errorTotalAirspeed2(i1:i2).'],'Edges',{hWbins,hWbins});
figure;
pcolor(hWbins,hWbins,hwerr.')
shading flat;
xlabel('Wind Difference [m/s]')
ylabel('TAS Difference [m/s]')
title(ncfilenameB,'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_hist_Diff_TAS_Wind', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_hist_Diff_TAS_Wind', string(timeWL(i1)),'_to_', string(timeWL(i2)),'sec', FileTag, '.png'));
end

hWbins = linspace(-3,3,200);
hwzerr = hist3([(wRD(3,i1:i2)-wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',errorTotalAirspeed2(i1:i2).'],'Edges',{hWbins,hWbins});
% hwzerr = hist3([wRD(3,i1:i2).',errorTotalAirspeed2(i1:i2).'],'Edges',{hWbins,hWbins});
figure;
pcolor(hWbins,hWbins,hwzerr.')
shading flat;
xlabel('Vertical Wind Difference [m/s]')
ylabel('TAS Difference [m/s]')
title(ncfilenameB,'Interpreter','none');

hWbins = linspace(-3,3,200);
% hwzerr = hist3([(wRD(3,i1:i2)-wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).',errorTotalAirspeed2(i1:i2).'],'Edges',{hWbins,hWbins});
hwzerr = hist3([wRD(3,i1:i2).',(wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1)).'],'Edges',{hWbins,hWbins});
figure;
pcolor(hWbins,hWbins,hwzerr.')
shading flat;
hold on;
plot(hWbins,hWbins,'g--')
xlabel('Radome Vertical Wind [m/s]')
ylabel('LAMS Vertical Wind [m/s]')
title(ncfilenameB,'Interpreter','none');


figure;
scatter3(wRD(3,i1:i2),wLAMS(3,i1:i2).*prod(foundNan(:,i1:i2),1),errorTotalAirspeed2(i1:i2),1,errorTotalAirspeed2(i1:i2));

figure;
scatter(((sqrt(sum(wRD(1:2,i1:i2).^2,1))-sqrt(sum(wLAMS(1:2,i1:i2).^2,1))).*prod(foundNan(:,i1:i2),1)),errorTotalAirspeed2(i1:i2))

figure;
xSize = 25;
ySize = 10;
xLeft = 0;
yTop = 3;
set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
set(gcf,'Position',[500 500 xSize*50 ySize*50]);
plot(timeWL/3600,sqrt(sum(wLAMS(1:2,:).^2,1)).*prod(foundNan,1),'.');
hold on;
plot(timeWL/3600,wLAMS(3,:).*prod(foundNan,1),'.');
hold on;
set(gca, 'ColorOrderIndex',1);
plot(timeWL/3600,[sqrt(vic(iD0-1+(1:DataLen)).^2+uic(iD0-1+(1:DataLen)).^2) wic(iD0-1+(1:DataLen))].','-')
xlabel('Time [h-UTC]')
ylabel('wind speed [m/s]')
title(ncfilenameB,'Interpreter','none');
grid on;
legend('LAMS Horizontal','LAMS Vertical','Radome Horizontal','Radome Vertical')
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_hv_LAMS_Winds', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_hv_LAMS_Winds', FileTag, '.png'));
end


figure;
xSize = 25;
ySize = 10;
xLeft = 0;
yTop = 3;
set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
set(gcf,'Position',[500 500 xSize*50 ySize*50]);
plot(timeWL/3600,(wLAMS.*prod(foundNan,1)).*[-1;-1;1],'.');
hold on;
set(gca, 'ColorOrderIndex',1);
plot(timeWL/3600,[vic(iD0-1+(1:DataLen)) uic(iD0-1+(1:DataLen)) wic(iD0-1+(1:DataLen))].','-')
xlabel('Time [h-UTC]')
ylabel('wind speed [m/s]')
title(ncfilenameB,'Interpreter','none');
grid on;
legend('LAMS VIC','LAMS UIC','LAMS WIC','VIC','UIC','WIC')
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Winds', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Winds', FileTag, '.png'));
end


% Plot the running average of fraction of peaks found for each beam
% The calculation is performed over 101 data points
DetFrac = zeros(size(foundPk));
DetFrac(1,:) = smooth(foundPk(1,:),101);  % Compute the average peak detection percentages
DetFrac(2,:) = smooth(foundPk(2,:),101);  % Compute the average peak detection percentages
DetFrac(3,:) = smooth(foundPk(3,:),101);  % Compute the average peak detection percentages
DetFracTot = smooth(prod(foundPk,1),101);

% DetStd = smooth_std(aDiff.*foundNan,101);  % Compute the average peak detection percentages
figure;
plot(timeWL/3600,DetFrac.','linewidth',1.5);
hold on;
plot(timeWL/3600,DetFracTot,'k--','linewidth',1.5);
grid on;
xlabel('Time [h-UTC]');
ylabel('Fraction of Peaks Detected');
legend('Beam 1','Beam 2','Beam 3','Total')
title(ncfilenameB,'Interpreter','none');
if SaveFigs
   fnum = gcf;
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Detected_Fraction', FileTag, '.fig'));
   saveas(fnum,strcat(SaveDirectory, 'Plots/', ncfilenameB(1:end-3), '_LAMS_Detected_Fraction', FileTag, '.png'));
end

if ~isempty(UHSAS_Suffix)
    figure; 
    scatter(uhsasConc(iD0-1+(1:DataLen)),uhsasDbar(iD0-1+(1:DataLen)),1,DetFracTot)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    grid on;
    xlabel('UHSAS Particle Concentration [cm^{3}]')
    ylabel('UHSAS Mean Particle Diameter [\mum]')
    colorbar;
    
    figure; 
    scatter3(uhsasConc(iD0-1+(1:DataLen)),uhsasDbar(iD0-1+(1:DataLen)),galt(iD0-1+(1:DataLen)),1,DetFracTot)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    grid on;
    xlabel('UHSAS Particle Concentration [cm^{3}]')
    ylabel('UHSAS Mean Particle Diameter [\mum]')
    zlabel('Altitude [m]')
    colorbar;
end