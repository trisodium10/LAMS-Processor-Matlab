

% ncfilenameB = 'ARISTO2017_rf01.nc';  % LAMS Data file
% ncfilenameO = 'ARISTO2017_rf01.nc';  % Aircraft Data File
% UHSAS_Suffix = '';  % Suffix for UHSAS variable names if it is in the instrument payload
% 
% Aircraft = 'GV';  % Use if this is a GV flight
% % Aircraft = 'C-130'; % Use if this is a C-130 flight
% 
% FileTag = '';
% 
% addpath('/scr/raf_data/ARISTO2017');

LoadData = 1;
LAMStimeshift = 0;  % used to adjust time latency between LAMS and aircraft data system
INStimeshift = 0;  % used to adjust time latency between Honeywell and aircraft data system
CMtimeshift = 0;  % used to adjust time latency between SDN-500 and aircraft data system

croll_offset = 0;  % fixed offset in roll for SDN-500

SubSet = 1; % only load a portion of the flight based on istart and istop
SetTime = 0; % if Subset =1, SetTime=1 uses time data to select the subset.  Setime=0, uses TAS to find the start and endpoints
TAS_TH = 50;  % TAS threshold for deciding when to start and stop subset selection in m/s

if strcmp(Aircraft,'GV')
    BeamList = [1; 1; 1; 0];  % List indicating what beams are used (1 for used, 0 for not used)
else
    BeamList = [1; 1; 1; 1];  % List indicating what beams are used (1 for used, 0 for not used)
end



load('BeamPointing.mat')





% SDN-500 location relative to Honeywell (Aircraft coord)
if strcmp(Aircraft,'GV')
    rL = [-10.305; -6.319; -1.359];  % GV
else
    rL = [416; 1443; 101]/100;  % C-130
end



if LoadData
    
    time = NCloaddata('Time',ncfilenameB,nan);
    
    time1 = time;
    
    time50 = time(1)+(1:(length(time)*50))*mean(diff(time))/50;
    time50 = time50(:);
    % time = time50;  % use for high rate data
    
    
    
    % Load LAMS Spectra and determine data rate in use
    if BeamList(1)
        lams_spectra1 = NCloadhist('BEAM1_LAMS',ncfilenameB);
        SpecSize = size(lams_spectra1);
        if SpecSize(2) == length(time50)
            time = time50;
            HRdata = 1;
            FileTag = [FileTag '_HR'];
        else
            FileTag = [FileTag '_1Hz'];
            HRdata = 0;
        end
        time = time+LAMStimeshift;
        Tlen = length(time);

        if SubSet
            if SetTime
                [~,istart] = min(abs(time-tstart));
                [~,istop] = min(abs(time-tstop));
            else
                tas = NCloaddata('TASX',ncfilenameO,time);
                istart = find(tas > TAS_TH,1);
                istop = find(tas > TAS_TH,1,'last');
            end
        end


        if SubSet == 0
            istart = 1;
            istop = Tlen;
            tstart = time(1);
            tstop = time(end);
        end

        if istop > Tlen
           istop = Tlen;
        end

        if SubSet
            time = time(istart:istop);
        end

        if SubSet
            lams_spectra1 = lams_spectra1(:,istart:istop);
        end
    end
    
    if BeamList(2)
        lams_spectra2 = NCloadhist('BEAM2_LAMS',ncfilenameB);
        if SubSet
            lams_spectra2 = lams_spectra2(:,istart:istop);
        end
    end
    
    if BeamList(3)
        lams_spectra3 = NCloadhist('BEAM3_LAMS',ncfilenameB);
        if SubSet
            lams_spectra3 = lams_spectra3(:,istart:istop);
        end
    end
    
    if BeamList(4)
        lams_spectra4 = NCloadhist('BEAM4_LAMS',ncfilenameB);
        if SubSet
            lams_spectra4 = lams_spectra4(:,istart:istop);
        end
    end
    
    % netcdf.close(ncid);
    
    % varid_speed = netcdf.inqVarID(ncid,'BEAM1_speed');
    % aSpeed1 = double(netcdf.getVar(ncid,varid_speed));
    % aSpeed1 = aSpeed1(istart:istop);
    %
    % varid_speed = netcdf.inqVarID(ncid,'BEAM2_speed');
    % aSpeed2 = double(netcdf.getVar(ncid,varid_speed));
    % aSpeed2 = aSpeed2(istart:istop);
    %
    % varid_speed = netcdf.inqVarID(ncid,'BEAM3_speed');
    % aSpeed3 = double(netcdf.getVar(ncid,varid_speed));
    % aSpeed3 = aSpeed3(istart:istop);
    
    % estimate time lag in LAMS data
    gglag = NCloaddata('GGRepLag_UDP',ncfilenameO,time);
    
%     timeS = time - gglag;
    timeS = time;
    
    tas = NCloaddata('TASX',ncfilenameO,timeS);
    thdg = NCloaddata('THDG',ncfilenameO,timeS+INStimeshift);
    pitch = NCloaddata('PITCH',ncfilenameO,timeS+INStimeshift);
    roll = NCloaddata('ROLL',ncfilenameO,timeS+INStimeshift);
    cthdg = NCloaddata('CTHDG_LAMS',ncfilenameB,timeS-CMtimeshift);
    cpitch = NCloaddata('CPITCH_LAMS',ncfilenameB,timeS-CMtimeshift);
    croll = NCloaddata('CROLL_LAMS',ncfilenameB,timeS-CMtimeshift)+croll_offset;
    calt = NCloaddata('CALT_LAMS',ncfilenameB,timeS-CMtimeshift);
    clon = NCloaddata('CLON_LAMS',ncfilenameO,timeS-CMtimeshift);
    clat = NCloaddata('CLAT_LAMS',ncfilenameO,timeS-CMtimeshift);
    cfominf = NCloaddata('FOMINF_LAMS',ncfilenameO,timeS-CMtimeshift);  % Figure of merrit for SDN-500, lower is better
    cmode = NCloaddata('MODE_LAMS',ncfilenameO,timeS-CMtimeshift);  % SDN-500 mode - 7 is normal operation
    cexvelerr = NCloaddata('EXVELERR_LAMS',ncfilenameO,timeS-CMtimeshift);  % SDN-500 expected velocity error 
    
    wic = NCloaddata('WIC',ncfilenameO,timeS);
    uic = NCloaddata('UIC',ncfilenameO,timeS);
    vic = NCloaddata('VIC',ncfilenameO,timeS);
    wsc = NCloaddata('WSC',ncfilenameO,timeS);
    wdc = NCloaddata('WDC',ncfilenameO,timeS);
    vnsc = NCloaddata('GGVNS',ncfilenameB,timeS);
    vewc = NCloaddata('GGVEW',ncfilenameB,timeS);
    attack = NCloaddata('ATTACK',ncfilenameO,timeS);
    attackrd = NCloaddata('AKRD',ncfilenameO,timeS);
    ssrd = NCloaddata('SSRD',ncfilenameO,timeS);
    gspd = NCloaddata('GGSPD',ncfilenameB,timeS);
    galt = NCloaddata('GGALT',ncfilenameB,timeS);
    gvspd = NCloaddata('GGVSPD',ncfilenameB,timeS);
    lat = NCloaddata('GGLAT',ncfilenameB,timeS);
    lon = NCloaddata('GGLON',ncfilenameB,timeS);
    
    adifr = NCloaddata('ADIFR',ncfilenameB,timeS);
    bdifr = NCloaddata('BDIFR',ncfilenameB,timeS);
    qcf = NCloaddata('QCF',ncfilenameB,timeS);
    mach = NCloaddata('MACHX',ncfilenameB,timeS);
    if strcmp(Aircraft,'C-130')
        psf = NCloaddata('PSFD',ncfilenameB,timeS);    %C-130
    else
        psf = NCloaddata('PSF',ncfilenameB,timeS);  % GV
    end
    ewx = NCloaddata('EWX',ncfilenameB,timeS);
    atx = NCloaddata('ATX',ncfilenameB,timeS);  % air temperature
    
    cvspd = -1.0*NCloaddata('CVSPD_LAMS',ncfilenameB,timeS);
    cvew = NCloaddata('CVEW_LAMS',ncfilenameB,timeS);
    cvns = NCloaddata('CVNS_LAMS',ncfilenameB,timeS);
    
%     UHSAS params when available
    if ~isempty(UHSAS_Suffix)
        uhsasConc500 = NCloaddata(['CONCU500' UHSAS_Suffix],ncfilenameO,timeS);
        uhsasConc100 = NCloaddata(['CONCU100' UHSAS_Suffix],ncfilenameO,timeS);
        uhsasConc = NCloaddata(['CONCU' UHSAS_Suffix],ncfilenameO,timeS);
        uhsasDbar = NCloaddata(['DBARU' UHSAS_Suffix],ncfilenameO,timeS);
    end
    
%     conc100 = NCloaddata('CONC1DC100_RWO',ncfilenameB,timeS);  % CDP concentration
    % tcab = NCloaddata('TCAB',ncfilenameO,time);
    
end


% Manual Pod pointing offsets (between beams and SDN-500)
if strcmp(Aircraft,'GV')
    RollPod = 0;  % pod roll depends on aircraft.  GV (0) and C-130 (pi)
    PitchPod = -0.0015; %-0.0026/2;% -0.2352*pi/180/2 % -0.2152*pi/180
    YawPod = 0;%0.0018;  % 0.2*pi/180
else
    RollPod = pi;
    PitchPod = 0; %0.1752*pi/180/2; %-0.2352*pi/180/2; % -0.2152*pi/180
end
Beam2ThetaAdj = 0*pi/180;
rB2v = cross(vC(:,2),[1;0;0]);
rB2v = rB2v/norm(rB2v);
vC2a = Rprt(Beam2ThetaAdj,rB2v)*vC(:,2);
k1L = RPYmat(RollPod,PitchPod,YawPod)*vC(:,1);
k4L = RPYmat(RollPod,PitchPod,YawPod)*vC(:,2);     % beam 4 points under the nose
k2L = RPYmat(RollPod,PitchPod,YawPod)*vC(:,4);  % beam 2 is forward pointing
k3L = RPYmat(RollPod,PitchPod,YawPod)*vC(:,3);
k1L(2:3) = -k1L(2:3);
k2L(2:3) = -k2L(2:3);
k3L(2:3) = -k3L(2:3);
k4L(2:3) = -k4L(2:3);

