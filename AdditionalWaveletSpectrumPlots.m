
% Define time limts of desired spectrum
% [~,i1] = min(abs(timeWL/3600-22.78));
% [~,i2] = min(abs(timeWL/3600-22.98));
% i2 = i2-2;

% [~,i1] = min(abs(timeWL/3600-21.91));
% [~,i2] = min(abs(timeWL/3600-21.94));

% [~,i1] = min(abs(timeWL/3600-22.5));
% [~,i2] = min(abs(timeWL/3600-22.7));
% i2 = i2-2;

% [~,i1] = min(abs(timeWL/3600-22.43));
% [~,i2] = min(abs(timeWL/3600-22.44));

% [~,i1] = min(abs(timeWL/3600-16.5));
% [~,i2] = min(abs(timeWL/3600-16.55));

% i1 = 1;
% % i2 = length(timeWL);
% [~,i2] = min(abs(timeWL/3600-22.8));

% [~,i1] = min(abs(timeWL/3600-22.16));
% [~,i2] = min(abs(timeWL/3600-22.18));

% i1 = 900;
% i2 = 1500;

SaveFigure = 0;

t_indices = (i1:i2);
t_indices2 = (i1:(i2-1));

NumBeams = sum(BeamList);

% Use to turn on different plotting options
PeakLocations = 1;  % Plot locations of all peaks found in the spectrum
PeakLocations_OnlyCurrent = 0;  % if PeakLocations==1, Only plot all the found peaks for the current instance
RadomeTrack = 1;   % Plot the location of the expected radome track
AcceptedPeaks = 1;  % Plot the peaks that are accepted and assumed
AllAcceptedPeaks = 1;  % If AcceptedPeaks==1, Plot the accepted/assumed peaks of the latest spectrum

% Plot Dimension adjustments
SpdWid = 40;  % Set the y Limit width
xSpace = 0;  % Extra Space left on x axis - often set to 0.35 to leave room for the legend

% Define markers to be used for the plot overlays
Mark_PeakLocation = 'bx';
Mark_AcceptedPeak = 'gx';
Mark_AssumedPeak = 'rx';
Mark_NewPeak = 'ko';
Mark_Radome = 'k:';



% figure size parameters
xSize = 25;
ySize = 10;
xLeft = 0;
yTop = 3;

LegendList = {};
LegendFigs = [];
LegendIndex = 1;

figure;
set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
set(gcf,'Position',[500 500 xSize*50 ySize*50]);
subplot(NumBeams,1,1);

lams_speed1 = [WLspec1(:,t_indices);WLspec1(end:-1:1,t_indices);WLspec1(:,t_indices);WLspec1(end:-1:1,t_indices)];
spdvect = ((1:size(lams_speed1,1))-1)*dfs*lambda/2;
% indexvect = (1:size(lams_speed1,2))+iD0+istart;
indexvect = timeWL(t_indices).'/3600;
indexvect2 = timeWL(t_indices2).'/3600;
h1 = pcolor(indexvect,spdvect,real(log10(lams_speed1))); 
LegendList{LegendIndex} = 'Spectrum';
LegendFigs(LegendIndex) = h1;
LegendIndex = LegendIndex +1;
shading flat;
hold on; 
% title(['Beam 1 ' ncfilenameB])
% plot(time/3600,iSpec(1,:),'k-');
% plot(time/3600,iTAS(1,(:,iD0+DataLen-1)),'k:');
% xlabel('Time [h UT]');
% ylabel('Array Index');
if PeakLocations
    if PeakLocations_OnlyCurrent
        h2 = plot(indexvect(end)*ones(size(spdSet{1,t_indices(end)},1)-1,1),spdSet{1,t_indices(end)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
    else
        for x1 = 1:length(t_indices)
            h2a = plot(indexvect(x1)*ones(size(spdSet{1,t_indices(x1)},1)-1,1),spdSet{1,t_indices(x1)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
            if ~isempty(h2a)
                h2 = h2a;
            end
        %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
        end
    end
    LegendList{LegendIndex} = 'Peaks';
    LegendFigs(LegendIndex) = h2;
    LegendIndex = LegendIndex +1;
end
if AcceptedPeaks
    if AllAcceptedPeaks
        aLAMSValid = aLAMSlist(1,t_indices);
        aLAMSValid(foundPk(1,t_indices)==0) = nan;
        aLAMSNotValid = aLAMSlist(1,t_indices);
        aLAMSNotValid(foundPk(1,t_indices)==1) = nan;
        h3 = plot(indexvect,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    else
        aLAMSValid = aLAMSlist(1,t_indices2);
        aLAMSValid(foundPk(1,t_indices2)==0) = nan;
        aLAMSNotValid = aLAMSlist(1,t_indices2);
        aLAMSNotValid(foundPk(1,t_indices2)==1) = nan;
        h3 = plot(indexvect2,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect2,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    end
end
if RadomeTrack
    h5 = plot(indexvect2,aLAMSlist(1,t_indices2)-aDiff(1,t_indices2),Mark_Radome,'linewidth',1.5);
    h6 = plot(indexvect(end),aLAMSlist(1,t_indices(end))-aDiff(1,t_indices(end)),Mark_NewPeak,'linewidth',1.5);
    LegendList{LegendIndex} = 'Estimated Radome Track';
    LegendFigs(LegendIndex) = h5;
    LegendIndex = LegendIndex +1;
    LegendList{LegendIndex} = 'Estimated New Peak Location';
    LegendFigs(LegendIndex) = h6;
    LegendIndex = LegendIndex +1;
end
% title(['Beam 1 ' ncfilenameB],'Interpreter', 'none')
title('Beam 1');
ylabel('LOS speed [m/s]');
% ylim([150,170]);
yMid1 = nanmean(aLAMSlist(1,t_indices2)-aDiff(1,t_indices2));
ylim([yMid1-SpdWid/2,yMid1+SpdWid/2]);
xlim([xlim*[1;0],xlim*[-1;1]*xSpace+xlim*[0;1]]);
% xlabel('1 Hz index');
% xlabel('Time [h-UTC]');
legend(LegendFigs,LegendList);



subplot(NumBeams,1,2);
LegendList = {};
LegendFigs = [];
LegendIndex = 1;

lams_speed2 = [WLspec2(:,t_indices);WLspec2(end:-1:1,t_indices);WLspec2(:,t_indices);WLspec2(end:-1:1,t_indices)];
spdvect = ((1:size(lams_speed1,1))-1)*dfs*lambda/2;
% indexvect = (1:size(lams_speed1,2))+iD0+istart;
indexvect = timeWL(t_indices).'/3600;
indexvect2 = timeWL(t_indices2).'/3600;
h1 = pcolor(indexvect,spdvect,real(log10(lams_speed2))); 
LegendList{LegendIndex} = 'Spectrum';
LegendFigs(LegendIndex) = h1;
LegendIndex = LegendIndex +1;
shading flat;
hold on; 
if PeakLocations
    if PeakLocations_OnlyCurrent
        h2 = plot(indexvect(end)*ones(size(spdSet{2,t_indices(end)},1)-1,1),spdSet{2,t_indices(end)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
    else
        for x1 = 1:length(t_indices)
            h2a = plot(indexvect(x1)*ones(size(spdSet{2,t_indices(x1)},1)-1,1),spdSet{2,t_indices(x1)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
            if ~isempty(h2a)
                h2 = h2a;
            end
        %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
        end
    end
    LegendList{LegendIndex} = 'Peaks';
    LegendFigs(LegendIndex) = h2;
    LegendIndex = LegendIndex +1;
end
if AcceptedPeaks
    if AllAcceptedPeaks
        aLAMSValid = aLAMSlist(2,t_indices);
        aLAMSValid(foundPk(2,t_indices)==0) = nan;
        aLAMSNotValid = aLAMSlist(2,t_indices);
        aLAMSNotValid(foundPk(2,t_indices)==1) = nan;
        h3 = plot(indexvect,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    else
        aLAMSValid = aLAMSlist(2,t_indices2);
        aLAMSValid(foundPk(2,t_indices2)==0) = nan;
        aLAMSNotValid = aLAMSlist(2,t_indices2);
        aLAMSNotValid(foundPk(2,t_indices2)==1) = nan;
        h3 = plot(indexvect2,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect2,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    end
end
if RadomeTrack
    h5 = plot(indexvect2,aLAMSlist(2,t_indices2)-aDiff(2,t_indices2),Mark_Radome,'linewidth',1.5);
    h6 = plot(indexvect(end),aLAMSlist(2,t_indices(end))-aDiff(2,t_indices(end)),Mark_NewPeak,'linewidth',1.5);
    LegendList{LegendIndex} = 'Estimated Radome Track';
    LegendFigs(LegendIndex) = h5;
    LegendIndex = LegendIndex +1;
    LegendList{LegendIndex} = 'Estimated New Peak Location';
    LegendFigs(LegendIndex) = h6;
    LegendIndex = LegendIndex +1;
end
title('Beam 2');
ylabel('LOS speed [m/s]');
% ylim([180,200]);
% ylim([120,180]);
yMid2 = nanmean(aLAMSlist(2,t_indices2)-aDiff(2,t_indices2));
ylim([yMid2-SpdWid/2,yMid2+SpdWid/2]);
xlim([xlim*[1;0],xlim*[-1;1]*xSpace+xlim*[0;1]]);
% xlabel('1 Hz index');
% xlabel('Time [h-UTC]');
% legend(LegendFigs,LegendList);

subplot(NumBeams,1,3);
LegendList = {};
LegendFigs = [];
LegendIndex = 1;

lams_speed3 = [WLspec3(:,t_indices);WLspec3(end:-1:1,t_indices);WLspec3(:,t_indices);WLspec3(end:-1:1,t_indices)];
spdvect = ((1:size(lams_speed3,1))-1)*dfs*lambda/2;
% indexvect = (1:size(lams_speed1,2))+iD0+istart;
indexvect = timeWL(t_indices).'/3600;
indexvect2 = timeWL(t_indices2).'/3600;
h1 = pcolor(indexvect,spdvect,real(log10(lams_speed3))); 
LegendList{LegendIndex} = 'Spectrum';
LegendFigs(LegendIndex) = h1;
LegendIndex = LegendIndex +1;
shading flat;
hold on; 
if PeakLocations
    if PeakLocations_OnlyCurrent
        h2 = plot(indexvect(end)*ones(size(spdSet{3,t_indices(end)},1)-1,1),spdSet{3,t_indices(end)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
    else
        for x1 = 1:length(t_indices)
            h2a = plot(indexvect(x1)*ones(size(spdSet{3,t_indices(x1)},1)-1,1),spdSet{3,t_indices(x1)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
            if ~isempty(h2a)
                h2 = h2a;
            end
        %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
        end
    end
    LegendList{LegendIndex} = 'Peaks';
    LegendFigs(LegendIndex) = h2;
    LegendIndex = LegendIndex +1;
end
if AcceptedPeaks
    if AllAcceptedPeaks
        aLAMSValid = aLAMSlist(3,t_indices);
        aLAMSValid(foundPk(3,t_indices)==0) = nan;
        aLAMSNotValid = aLAMSlist(3,t_indices);
        aLAMSNotValid(foundPk(3,t_indices)==1) = nan;
        h3 = plot(indexvect,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    else
        aLAMSValid = aLAMSlist(3,t_indices2);
        aLAMSValid(foundPk(3,t_indices2)==0) = nan;
        aLAMSNotValid = aLAMSlist(3,t_indices2);
        aLAMSNotValid(foundPk(3,t_indices2)==1) = nan;
        h3 = plot(indexvect2,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
        h4 = plot(indexvect2,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
        LegendList{LegendIndex} = 'Accepted Peak';
        LegendFigs(LegendIndex) = h3;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Assumed Peak';
        LegendFigs(LegendIndex) = h4;
        LegendIndex = LegendIndex +1;  
    end
end
if RadomeTrack
    h5 = plot(indexvect2,aLAMSlist(3,t_indices2)-aDiff(3,t_indices2),Mark_Radome,'linewidth',1.5);
    h6 = plot(indexvect(end),aLAMSlist(3,t_indices(end))-aDiff(3,t_indices(end)),Mark_NewPeak,'linewidth',1.5);
    LegendList{LegendIndex} = 'Estimated Radome Track';
    LegendFigs(LegendIndex) = h5;
    LegendIndex = LegendIndex +1;
    LegendList{LegendIndex} = 'Estimated New Peak Location';
    LegendFigs(LegendIndex) = h6;
    LegendIndex = LegendIndex +1;
end
title('Beam 3');
ylabel('LOS speed [m/s]');
% ylim([150,170]);
% ylim([90,150]);
yMid3 = nanmean(aLAMSlist(3,t_indices2)-aDiff(3,t_indices2));
ylim([yMid3-SpdWid/2,yMid3+SpdWid/2]);
xlim([xlim*[1;0],xlim*[-1;1]*xSpace+xlim*[0;1]]);
if NumBeams == 3
    % xlabel('1 Hz index');
    xlabel('Time [h-UTC]');
end


if NumBeams == 4
   subplot(NumBeams,1,4);
    LegendList = {};
    LegendFigs = [];
    LegendIndex = 1;

    lams_speed4 = [WLspec4(:,t_indices);WLspec4(end:-1:1,t_indices);WLspec4(:,t_indices);WLspec4(end:-1:1,t_indices)];
    spdvect = ((1:size(lams_speed4,1))-1)*dfs*lambda/2;
    % indexvect = (1:size(lams_speed1,2))+iD0+istart;
    indexvect = timeWL(t_indices).'/3600;
    indexvect2 = timeWL(t_indices2).'/3600;
    h1 = pcolor(indexvect,spdvect,real(log10(lams_speed4))); 
    LegendList{LegendIndex} = 'Spectrum';
    LegendFigs(LegendIndex) = h1;
    LegendIndex = LegendIndex +1;
    shading flat;
    hold on; 
    if PeakLocations
        if PeakLocations_OnlyCurrent
            h2 = plot(indexvect(end)*ones(size(spdSet{4,t_indices(end)},1)-1,1),spdSet{4,t_indices(end)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
        else
            for x1 = 1:length(t_indices)
                h2a = plot(indexvect(x1)*ones(size(spdSet{4,t_indices(x1)},1)-1,1),spdSet{4,t_indices(x1)}(2:end),Mark_PeakLocation,'linewidth',1.5); 
                if ~isempty(h2a)
                    h2 = h2a;
                end
            %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
            end
        end
        LegendList{LegendIndex} = 'Peaks';
        LegendFigs(LegendIndex) = h2;
        LegendIndex = LegendIndex +1;
    end
    if AcceptedPeaks
        if AllAcceptedPeaks
            aLAMSValid = aLAMSlist(4,t_indices);
            aLAMSValid(foundPk(4,t_indices)==0) = nan;
            aLAMSNotValid = aLAMSlist(4,t_indices);
            aLAMSNotValid(foundPk(4,t_indices)==1) = nan;
            h3 = plot(indexvect,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
            h4 = plot(indexvect,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
            LegendList{LegendIndex} = 'Accepted Peak';
            LegendFigs(LegendIndex) = h3;
            LegendIndex = LegendIndex +1;
            LegendList{LegendIndex} = 'Assumed Peak';
            LegendFigs(LegendIndex) = h4;
            LegendIndex = LegendIndex +1;  
        else
            aLAMSValid = aLAMSlist(4,t_indices2);
            aLAMSValid(foundPk(4,t_indices2)==0) = nan;
            aLAMSNotValid = aLAMSlist(4,t_indices2);
            aLAMSNotValid(foundPk(4,t_indices2)==1) = nan;
            h3 = plot(indexvect2,aLAMSValid,Mark_AcceptedPeak,'linewidth',2);
            h4 = plot(indexvect2,aLAMSNotValid,Mark_AssumedPeak,'linewidth',2);
            LegendList{LegendIndex} = 'Accepted Peak';
            LegendFigs(LegendIndex) = h3;
            LegendIndex = LegendIndex +1;
            LegendList{LegendIndex} = 'Assumed Peak';
            LegendFigs(LegendIndex) = h4;
            LegendIndex = LegendIndex +1;  
        end
    end
    if RadomeTrack
        h5 = plot(indexvect2,aLAMSlist(4,t_indices2)-aDiff(4,t_indices2),Mark_Radome,'linewidth',1.5);
        h6 = plot(indexvect(end),aLAMSlist(4,t_indices(end))-aDiff(4,t_indices(end)),Mark_NewPeak,'linewidth',1.5);
        LegendList{LegendIndex} = 'Estimated Radome Track';
        LegendFigs(LegendIndex) = h5;
        LegendIndex = LegendIndex +1;
        LegendList{LegendIndex} = 'Estimated New Peak Location';
        LegendFigs(LegendIndex) = h6;
        LegendIndex = LegendIndex +1;
    end
    title('Beam 4');
    ylabel('LOS speed [m/s]');
%     ylim([150,170]);
    yMid4 = nanmean(aLAMSlist(4,t_indices2)-aDiff(4,t_indices2));
    ylim([yMid4-SpdWid/2,yMid4+SpdWid/2]);
    xlim([xlim*[1;0],xlim*[-1;1]*xSpace+xlim*[0;1]]);

    % xlabel('1 Hz index');
    xlabel('Time [h-UTC]');
end

% lams_speed2 = [WLspec2(:,t_indices);WLspec2(end:-1:1,t_indices);WLspec2(:,t_indices);WLspec2(end:-1:1,t_indices)];
% spdvect = ((1:size(lams_speed2,1))-1)*dfs*lambda/2;
% % indexvect = (1:size(lams_speed2,2))+iD0+istart;
% indexvect = timeWL(t_indices).'/3600;
% h1 = pcolor(indexvect,spdvect,real(log10(lams_speed2))); 
% shading flat;
% hold on; 
% set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
% set(gcf,'Position',[500 500 xSize*50 ySize*50]);
% for x1 = 1:length(t_indices)
%     h2 = plot(indexvect(x1)*ones(size(spdSet{2,t_indices(x1)},1)-1,1),spdSet{2,t_indices(x1)}(2:end),'b.'); 
% %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
% end
% aLAMSValid = aLAMSlist(2,t_indices);
% aLAMSValid(foundPk(2,t_indices)==0) = nan;
% aLAMSNotValid = aLAMSlist(2,t_indices);
% aLAMSNotValid(foundPk(2,t_indices)==1) = nan;
% h3 = plot(indexvect,aLAMSValid,'g.','linewidth',2);
% h4 = plot(indexvect,aLAMSNotValid,'r.','linewidth',2);
% % plot(indexvect,aLAMSlist(2,:),'k--','linewidth',2);
% h5 = plot(indexvect,aLAMSlist(2,t_indices)-aDiff(2,t_indices),'k:','linewidth',1.5);
% title(['Beam 2 ' ncfilenameB],'Interpreter','none')
% ylabel('LOS speed [m/s]');
% % xlabel('1 Hz index');
% xlabel('Time [h-UTC]');
% legend([h1,h2,h3,h4,h5],{'Spectrum','Peaks','Accepted Peak','Assumed Peak','Estimated Radome Track'});
% 
% lams_speed3 = [WLspec3(:,t_indices);WLspec3(end:-1:1,t_indices);WLspec3(:,t_indices);WLspec3(end:-1:1,t_indices)];
% % lams_speed3 = [lams_spectra3(:,t_indices);lams_spectra3(end:-1:1,t_indices);lams_spectra3(:,t_indices);lams_spectra3(end:-1:1,t_indices)];
% spdvect = ((1:size(lams_speed3,1))-1)*dfs*lambda/2;
% % indexvect = (1:size(lams_speed3,2))+iD0+istart;
% indexvect = timeWL(t_indices).'/3600;
% figure; 
% h1 = pcolor(indexvect,spdvect,real(log10(lams_speed3))); 
% shading flat;
% ylim([130, 190]);
% xlim([22.78, 23]);
% hold on; 
% set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
% set(gcf,'Position',[500 500 xSize*50 ySize*50]);
% for x1 = 1:length(t_indices)
%     h2 = plot(indexvect(x1)*ones(size(spdSet{3,t_indices(x1)},1)-1,1),spdSet{3,t_indices(x1)}(2:end),'b.'); 
% %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
% end
% aLAMSValid = aLAMSlist(3,t_indices);
% aLAMSValid(foundPk(3,t_indices)==0) = nan;
% aLAMSNotValid = aLAMSlist(3,t_indices);
% aLAMSNotValid(foundPk(3,t_indices)==1) = nan;
% h3 = plot(indexvect,aLAMSValid,'g.','linewidth',2);
% h4 = plot(indexvect,aLAMSNotValid,'r.','linewidth',2);
% % plot(indexvect,aLAMSlist(3,:),'k--','linewidth',2);
% h5 = plot(indexvect,aLAMSlist(3,t_indices)-aDiff(3,t_indices),'k:','linewidth',1.5);
% title(['Beam 3 ' ncfilenameB],'Interpreter','none')
% ylabel('LOS speed [m/s]');
% % xlabel('1 Hz index');
% xlabel('Time [h-UTC]');
% legend([h1,h2,h3,h4,h5],{'Spectrum','Peaks','Accepted Peak','Assumed Peak','Estimated Radome Track'});
% 
% if BeamList(4)
%     lams_speed4 = [WLspec4(:,t_indices);WLspec4(end:-1:1,t_indices);WLspec4(:,t_indices);WLspec4(end:-1:1,t_indices)];
%     spdvect = ((1:size(lams_speed4,1))-1)*dfs*lambda/2;
%     % indexvect = (1:size(lams_speed4,2))+iD0+istart;
%     indexvect = timeWL(t_indices).'/3600;
%     figure; 
%     h1 = pcolor(indexvect,spdvect,real(log10(lams_speed3))); 
%     shading flat;
%     hold on; 
%     set(gcf,'PaperPosition',[xLeft, yTop, xSize, ySize]);
%     set(gcf,'Position',[500 500 xSize*50 ySize*50]);
%     for x1 = 1:length(t_indices)
%         h2 = plot(indexvect(x1)*ones(size(spdSet{4,t_indices(x1)},1)-1,1),spdSet{4,t_indices(x1)}(2:end),'b.'); 
%     %    plot(time(x1)/3600*ones(size(i1{x1})),i1{x1},'r.'); 
%     end
%     aLAMSValid = aLAMSlist(4,t_indices);
%     aLAMSValid(foundPk(4,t_indices)==0) = nan;
%     aLAMSNotValid = aLAMSlist(4,t_indices);
%     aLAMSNotValid(foundPk(4,t_indices)==1) = nan;
%     h3 = plot(indexvect,aLAMSValid,'g.','linewidth',2);
%     h4 = plot(indexvect,aLAMSNotValid,'r.','linewidth',2);
%     % plot(indexvect,aLAMSlist(3,:),'k--','linewidth',2);
%     h5 = plot(indexvect,aLAMSlist(4,t_indices)-aDiff(4,t_indices),'k:','linewidth',1.5);
%     title(['Beam 4 ' ncfilenameB],'Interpreter','none')
%     ylabel('LOS speed [m/s]');
%     % xlabel('1 Hz index');
%     xlabel('Time [h-UTC]');
%     legend([h1,h2,h3,h4,h5],{'Spectrum','Peaks','Accepted Peak','Assumed Peak','Estimated Radome Track'});
% end


if SaveFigure
    saveas(gcf,['SpectrumPeakTracks_' ncfilenameB(1:end-3) '_T' num2str(round(timeWL(i1))) 'sUTC_T' num2str(round(timeWL(i2))) 'sUTC.fig']);
    saveas(gcf,['SpectrumPeakTracks_' ncfilenameB(1:end-3) '_T' num2str(round(timeWL(i1))) 'sUTC_T' num2str(round(timeWL(i2))) 'sUTC.png']);
end