

% % Fill in full nc file time series
% t_nc = time1;
% tsize = size(t_nc);
% tLen = length(t_nc);
% iOffset = iD0+istart-2;

% Fill in only Wavelet processed data
t_nc = timeWL;
tsize = size(t_nc);
tLen = length(t_nc);
iOffset = 0;



attack_L = zeros(tsize);
attack_L(iOffset+(1:DataLen)) = attackLAMS*180/pi;

ss_L = zeros(tsize);
ss_L(iOffset+(1:DataLen)) = ssLAMS*180/pi;

tas_L = zeros(tsize);
tas_L(iOffset+(1:DataLen)) = tasLAMS;

tas_L_U = zeros(tsize);
tas_L_U(iOffset+(1:DataLen)) = sqrt(var_tasLAMS);

if BeamList(1)
    v_LOS1 = zeros(tsize);
    v_LOS1(iOffset+(1:DataLen)) = aLAMSlist(1,:);

    fPk1 = zeros(tsize);
    fPk1(iOffset+(1:DataLen)) = foundPk(1,:);

    v_LOS1_U = zeros(tsize);
    v_LOS1_U(iOffset+(1:DataLen)) = sqrt(covLAMSLOSStat(1,:));
end

if BeamList(2)
    v_LOS2 = zeros(tsize);
    v_LOS2(iOffset+(1:DataLen)) = aLAMSlist(2,:);

    fPk2 = zeros(tsize);
    fPk2(iOffset+(1:DataLen)) = foundPk(2,:);

    v_LOS2_U = zeros(tsize);
    v_LOS2_U(iOffset+(1:DataLen)) = sqrt(covLAMSLOSStat(2,:));
end

if BeamList(3)
    v_LOS3 = zeros(tsize);
    v_LOS3(iOffset+(1:DataLen)) = aLAMSlist(3,:);

    fPk3 = zeros(tsize);
    fPk3(iOffset+(1:DataLen)) = foundPk(3,:);

    v_LOS3_U = zeros(tsize);
    v_LOS3_U(iOffset+(1:DataLen)) = sqrt(covLAMSLOSStat(3,:));
end

if BeamList(4)
    v_LOS4 = zeros(tsize);
    v_LOS4(iOffset+(1:DataLen)) = aLAMSlist(4,:);

    fPk4 = zeros(tsize);
    fPk4(iOffset+(1:DataLen)) = foundPk(4,:);

    v_LOS4_U = zeros(tsize);
    v_LOS4_U(iOffset+(1:DataLen)) = sqrt(covLAMSLOSStat(4,:));
end
v_L_AC_X = zeros(tsize);
v_L_AC_X(iOffset+(1:DataLen)) = vLAMSac(1,:);

v_L_AC_Y = zeros(tsize);
v_L_AC_Y(iOffset+(1:DataLen)) = vLAMSac(2,:);

v_L_AC_Z = zeros(tsize);
v_L_AC_Z(iOffset+(1:DataLen)) = vLAMSac(3,:);

v_L_AC_X_U = zeros(tsize);
v_L_AC_X_U(iOffset+(1:DataLen)) = sqrt(covLAMSacINS(1,:));

v_L_AC_Y_U = zeros(tsize);
v_L_AC_Y_U(iOffset+(1:DataLen)) = sqrt(covLAMSacINS(2,:));

v_L_AC_Z_U = zeros(tsize);
v_L_AC_Z_U(iOffset+(1:DataLen)) = sqrt(covLAMSacINS(3,:));

v_K_AC_X = zeros(tsize);
v_K_AC_X(iOffset+(1:DataLen)) = vRDEst(1,:);

v_K_AC_Y = zeros(tsize);
v_K_AC_Y(iOffset+(1:DataLen)) = vRDEst(2,:);

v_K_AC_Z = zeros(tsize);
v_K_AC_Z(iOffset+(1:DataLen)) = vRDEst(3,:);

v_K_AC_X_U = zeros(tsize);
v_K_AC_X_U(iOffset+(1:DataLen)) = sqrt(covLAMSac(1,:));

v_K_AC_Y_U = zeros(tsize);
v_KAC_Y_U(iOffset+(1:DataLen)) = sqrt(covLAMSac(2,:));

v_K_AC_Z_U = zeros(tsize);
v_K_AC_Z_U(iOffset+(1:DataLen)) = sqrt(covLAMSac(3,:));


v_L_L_X = zeros(tsize);
v_L_L_X(iOffset+(1:DataLen)) = vLAMSwl(1,:);

v_L_L_Y = zeros(tsize);
v_L_L_Y(iOffset+(1:DataLen)) = vLAMSwl(2,:);

v_L_L_Z = zeros(tsize);
v_L_L_Z(iOffset+(1:DataLen)) = vLAMSwl(3,:);

v_L_L_X_U = zeros(tsize);
v_L_L_X_U(iOffset+(1:DataLen)) = sqrt(covLAMSwlStat(1,:));

v_L_L_Y_U = zeros(tsize);
v_L_L_Y_U(iOffset+(1:DataLen)) = sqrt(covLAMSwlStat(2,:));

v_L_L_Z_U = zeros(tsize);
v_L_L_Z_U(iOffset+(1:DataLen)) = sqrt(covLAMSwlStat(3,:));

% winds in a global coordinate frame
wnsLAMS = zeros(tsize);
wnsLAMS(iOffset+(1:DataLen)) = -wLAMS(1,:);

wewLAMS = zeros(tsize);
wewLAMS(iOffset+(1:DataLen)) = -wLAMS(2,:);

wicLAMS= zeros(tsize);
wicLAMS(iOffset+(1:DataLen)) = wLAMS(3,:);

% wind uncertainties in global coordinate frame, excluding Kalman
% difference term
wnsLAMS_U = zeros(tsize);
wnsLAMS_U(iOffset+(1:DataLen)) = sqrt(covLAMSglob(1,:)+covLAMSglobINS(1,:));

wewLAMS_U = zeros(tsize);
wewLAMS_U(iOffset+(1:DataLen)) = sqrt(covLAMSglob(2,:)+covLAMSglobINS(2,:));

wicLAMS_U = zeros(tsize);
wicLAMS_U(iOffset+(1:DataLen)) = sqrt(covLAMSglob(3,:)+covLAMSglobINS(3,:));

wscLAMS = sqrt(wnsLAMS.^2+wewLAMS.^2);
wdcLAMS = atan2(-wewLAMS,-wnsLAMS)*180/pi;


TimeInterval = ncreadatt(ncfilenameB,'/','TimeInterval');
Platform = ncreadatt(ncfilenameB,'/','Platform');
FlightDate = ncreadatt(ncfilenameB,'/','FlightDate');
FlightNumber = ncreadatt(ncfilenameB,'/','FlightNumber');

if WriteNetCDF
    current_dir = cd();
    cd(SaveDirectory)
    nc_filename = [SaveDirectory ncfilenameB(1:end-3) '_LAMS' FileTag '.nc'];
    nccreate(nc_filename,'Time','Dimensions',{'time',tLen},'Datatype','double');
    ncwrite(nc_filename,'Time',t_nc);
    nccreate(nc_filename,'ATTACK_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'ATTACK_LAMS','units','degrees');
    ncwrite(nc_filename,'ATTACK_LAMS',attack_L);
    nccreate(nc_filename,'SS_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'SS_LAMS','units','degrees');
    ncwrite(nc_filename,'SS_LAMS',ss_L);
    
    nccreate(nc_filename,'TAS_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'TAS_LAMS','units','m/s');
    ncwrite(nc_filename,'TAS_LAMS',tas_L);
    nccreate(nc_filename,'TAS_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'TAS_LAMS_Uncertainty','description','Estimated Uncertainty in TAS_LAMS.');
    ncwriteatt(nc_filename,'TAS_LAMS_Uncertainty','units','m/s');
    ncwrite(nc_filename,'TAS_LAMS_Uncertainty',tas_L_U);
    
    
    if BeamList(1) == 1
        nccreate(nc_filename,'V_LOS_Beam1','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam1','units','m/s');
        ncwriteatt(nc_filename,'V_LOS_Beam1','description','Air velocity along the 1st LAMS beam.');
        ncwrite(nc_filename,'V_LOS_Beam1',v_LOS1);
        nccreate(nc_filename,'Beam1_Found_Peak','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'Beam1_Found_Peak','description','1 = a signal was found, 0 = no signal found');
        ncwrite(nc_filename,'Beam1_Found_Peak',fPk1);
        nccreate(nc_filename,'V_LOS_Beam1_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam1_Uncertainty','units','m/s');
        ncwrite(nc_filename,'V_LOS_Beam1_Uncertainty',v_LOS1_U);
    end
    
    if BeamList(2) == 1
        nccreate(nc_filename,'V_LOS_Beam2','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam2','units','m/s');
        ncwriteatt(nc_filename,'V_LOS_Beam2','description','Air velocity along the 2nd LAMS beam.');
        ncwrite(nc_filename,'V_LOS_Beam2',v_LOS2);
        nccreate(nc_filename,'Beam2_Found_Peak','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'Beam2_Found_Peak','description','1 = a signal was found, 0 = no signal found');
        ncwrite(nc_filename,'Beam2_Found_Peak',fPk2);
        nccreate(nc_filename,'V_LOS_Beam2_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam2_Uncertainty','units','m/s');
        ncwrite(nc_filename,'V_LOS_Beam2_Uncertainty',v_LOS2_U);
    end
    
    if BeamList(3) == 1
        nccreate(nc_filename,'V_LOS_Beam3','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam3','units','m/s');
        ncwriteatt(nc_filename,'V_LOS_Beam3','description','Air velocity along the 3rd LAMS beam.');
        ncwrite(nc_filename,'V_LOS_Beam3',v_LOS3);
        nccreate(nc_filename,'Beam3_Found_Peak','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'Beam3_Found_Peak','description','1 = a signal was found, 0 = no signal found');
        ncwrite(nc_filename,'Beam3_Found_Peak',fPk3);
        nccreate(nc_filename,'V_LOS_Beam3_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam3_Uncertainty','units','m/s');
        ncwrite(nc_filename,'V_LOS_Beam3_Uncertainty',v_LOS3_U);
    end
    
    if BeamList(4) == 1
        nccreate(nc_filename,'V_LOS_Beam4','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam4','units','m/s');
        ncwriteatt(nc_filename,'V_LOS_Beam4','description','Air velocity along the 4th LAMS beam.');
        ncwrite(nc_filename,'V_LOS_Beam4',v_LOS4);
        nccreate(nc_filename,'Beam4_Found_Peak','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'Beam4_Found_Peak','description','1 = a signal was found, 0 = no signal found');
        ncwrite(nc_filename,'Beam4_Found_Peak',fPk4);
        nccreate(nc_filename,'V_LOS_Beam4_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
        ncwriteatt(nc_filename,'V_LOS_Beam4_Uncertainty','units','m/s');
        ncwrite(nc_filename,'V_LOS_Beam4_Uncertainty',v_LOS4_U);
    end
    
    % Write out direct 3-D Measurements from LAMS
    nccreate(nc_filename,'AIRSPEED_X_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS','description','Air velocity along the forward direction in the aircraft coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_X_LAMS',v_L_AC_X);
    nccreate(nc_filename,'AIRSPEED_Y_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS','description','Air velocity along the right hand direction in the aircraft coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Y_LAMS',v_L_AC_Y);
    nccreate(nc_filename,'AIRSPEED_Z_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS','description','Air velocity along the downward direction in the aircraft coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Z_LAMS',v_L_AC_Z);
    
    nccreate(nc_filename,'AIRSPEED_X_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS_Uncertainty','units','m/s');
    ncwrite(nc_filename,'AIRSPEED_X_LAMS_Uncertainty',v_L_AC_X_U);
    nccreate(nc_filename,'AIRSPEED_Y_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS_Uncertainty','units','m/s');
    ncwrite(nc_filename,'AIRSPEED_Y_LAMS_Uncertainty',v_L_AC_Y_U);
    nccreate(nc_filename,'AIRSPEED_Z_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS_Uncertainty','units','m/s');
    ncwrite(nc_filename,'AIRSPEED_Z_LAMS_Uncertainty',v_L_AC_Z_U);
    
    % Write out Kalamin Filtered 3-D Measurements
    nccreate(nc_filename,'AIRSPEED_X_EST','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_X_EST','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_X_EST','description','Air velocity along the forward direction in the aircraft coordinate frame.  Estiamted using LAMS and Kalman Filter.');
    ncwrite(nc_filename,'AIRSPEED_X_EST',v_K_AC_X);
    nccreate(nc_filename,'AIRSPEED_Y_EST','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Y_EST','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Y_EST','description','Air velocity along the right hand direction in the aircraft coordinate frame.  Estiamted using LAMS and Kalman Filter.');
    ncwrite(nc_filename,'AIRSPEED_Y_EST',v_K_AC_Y);
    nccreate(nc_filename,'AIRSPEED_Z_EST','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Z_EST','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Z_EST','description','Air velocity along the downward direction in the aircraft coordinate frame.  Estiamted using LAMS and Kalman Filter.');
    ncwrite(nc_filename,'AIRSPEED_Z_EST',v_K_AC_Z);
    
    
    nccreate(nc_filename,'AIRSPEED_X_LAMS_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS_LAMS','description','Air velocity along the forward direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_X_LAMS_LAMS',v_L_L_X);
    nccreate(nc_filename,'AIRSPEED_Y_LAMS_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS_LAMS','description','Air velocity along the right hand direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Y_LAMS_LAMS',v_L_L_Y);
    nccreate(nc_filename,'AIRSPEED_Z_LAMS_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS_LAMS','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS_LAMS','description','Air velocity along the downward direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Z_LAMS_LAMS',v_L_L_Z);
    
    nccreate(nc_filename,'AIRSPEED_X_LAMS_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_X_LAMS_LAMS_Uncertainty','description','Statistical uncertainty in the air velocity along the forward direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_X_LAMS_LAMS_Uncertainty',v_L_L_X_U);
    nccreate(nc_filename,'AIRSPEED_Y_LAMS_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Y_LAMS_LAMS_Uncertainty','description','Statistical uncertainty in the air velocity along the right hand direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Y_LAMS_LAMS_Uncertainty',v_L_L_Y_U);
    nccreate(nc_filename,'AIRSPEED_Z_LAMS_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'AIRSPEED_Z_LAMS_LAMS_Uncertainty','description','Statistical uncertainty in the air velocity along the downward direction in the LAMS coordinate frame.');
    ncwrite(nc_filename,'AIRSPEED_Z_LAMS_LAMS_Uncertainty',v_L_L_Z_U);
    
    nccreate(nc_filename,'WIC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'WIC_LAMS','units','m/s');
    ncwriteatt(nc_filename,'WIC_LAMS','description','wind velocity in the vertical direction');
    ncwrite(nc_filename,'WIC_LAMS',wicLAMS);
    nccreate(nc_filename,'UIC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'UIC_LAMS','units','m/s');
    ncwriteatt(nc_filename,'UIC_LAMS','description','wind velocity in the east/west direction');
    ncwrite(nc_filename,'UIC_LAMS',wewLAMS);
    nccreate(nc_filename,'VIC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'VIC_LAMS','units','m/s');
    ncwriteatt(nc_filename,'VIC_LAMS','description','wind velocity in the north/south direction');
    ncwrite(nc_filename,'VIC_LAMS',wnsLAMS);
    nccreate(nc_filename,'WSC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'WSC_LAMS','units','m/s');
    ncwriteatt(nc_filename,'WSC_LAMS','description','wind magnitude');
    ncwrite(nc_filename,'WSC_LAMS',wscLAMS);
    nccreate(nc_filename,'WDC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'WDC_LAMS','units','degrees');
    ncwriteatt(nc_filename,'WDC_LAMS','description','wind direction');
    ncwrite(nc_filename,'WDC_LAMS',wdcLAMS);

    nccreate(nc_filename,'WIC_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'WIC_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'WIC_LAMS_Uncertainty','description','Statistical and INS uncertainty in the wind velocity in the vertical direction.');
    ncwrite(nc_filename,'WIC_LAMS_Uncertainty',wicLAMS);
    nccreate(nc_filename,'UIC_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'UIC_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'UIC_LAMS_Uncertainty','description','Statistical and INS uncertainty in the wind velocity in the east/west direction.');
    ncwrite(nc_filename,'UIC_LAMS_Uncertainty',wewLAMS);
    nccreate(nc_filename,'VIC_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
    ncwriteatt(nc_filename,'VIC_LAMS_Uncertainty','units','m/s');
    ncwriteatt(nc_filename,'VIC_LAMS_Uncertainty','description','Statistical and INS uncertainty in the wind velocity in the north/south direction.');
    ncwrite(nc_filename,'VIC_LAMS_Uncertainty',wnsLAMS);
%     nccreate(nc_filename,'WSC_LAMS_Uncertainty','Dimensions',{'time',tLen},'Datatype','double');
%     ncwrite(nc_filename,'WSC_LAMS_Uncertainty',wscLAMS);
%     nccreate(nc_filename,'WDC_LAMS','Dimensions',{'time',tLen},'Datatype','double');
%     ncwrite(nc_filename,'WDC_LAMS',wdcLAMS);
    
    ncwriteatt(nc_filename,'/','creation_time',datestr(now));
    
    ncwriteatt(nc_filename,'/','FlightDate',FlightDate);
    ncwriteatt(nc_filename,'/','FlightNumber',FlightNumber);
    ncwriteatt(nc_filename,'/','Platform',Platform);
    ncwriteatt(nc_filename,'/','TimeInterval',TimeInterval);
    cd(current_dir);
end