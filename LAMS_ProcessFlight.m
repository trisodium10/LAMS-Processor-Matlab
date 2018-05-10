% % ARISTO RF06
% iD0 = 30400;%9345; 30400
% DataLen = 400; %400; 400


% % % % % % % % % % % Full Flight
% iD0 = 1;% 2000
% DataLen = size(lams_spectra1,2)-iD0;
% BeamList = [1 1 1 0];


% % % ARISTO 2017 RF04
% [~,iD0] = min(abs(time/3600-20.08));% 2000
% % [~,b] = min(abs(time/3600-24.5));
% % [~,b] = min(abs(time/3600-23));  % 22.84
% % DataLen = b-iD0; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% DataLen = size(lams_spectra1,2)-iD0;
% BeamList = [1 1 1 0];


% % % ARISTO 2017 RF01, high rate
% [~,iD0] = min(abs(time/3600-22.4));% 2000
% [~,b] = min(abs(time/3600-22.5));
% DataLen = size(lams_spectra1,2)-iD0;
% BeamList = [1 1 1 0];
% 

% 
% Tl2aTset = eye(3);
% % Tl2aTset = Tl2aT;

% % ARISTO 2016 RF06 INS Operational
% [~,iD0] = min(abs(time/3600-19.065));
% % [~,b] = min(abs(time/3600-23.1));
% DataLen = size(lams_spectra1,2)-iD0; %b-iD0; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% % % [~,iD0] = min(abs(time/3600-22.1));
% % % DataLen = 500; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% BeamList = [1 1 1 0];

% % ARISTO RF01
% iD0 = 9345; % 7964
% DataLen = 400;

% % % % % ARISTO 2017 RF01
% iD0 = 6560;% 2000
% % DataLen = size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% DataLen = 6702-iD0;
% BeamList = [1 1 1 0];


% % % ARISTO 2017 RF04
% [~,iD0] = min(abs(time/3600-20.08));% 2000
% % [~,b] = min(abs(time/3600-24.5));
% % [~,b] = min(abs(time/3600-23));  % 22.84
% % DataLen = b-iD0; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% DataLen = size(lams_spectra1,2)-iD0;
% BeamList = [1 1 1 0];

% WriteNetCDF = 0;
RotationCorrect = 1;  % subtract the aircraft rotation component from the aircraft frame velocity

aWL = 1:5;  % wavelet widths to be used in peak finding

stdNoPk = 10;  % LOS standard deviation (in m/s) to be used when no peak is found

NoPkProb = exp(-3^2);  % derating by assuming the correct is not detected (derating defined by number of standard diviations)

PkRepeatDerate = 0.2;  % Factor used to derate probability of previous peaks - nominally 0.5

ResFactor = 1/5.0;  % Sub grid resolution factor for speed uncertainty.  Testing of peak finding functions suggests 1/5.

% ARISTO 2017 RF03 - testing leg
% [~,iD0] = min(abs(time/3600-22.1));
% [~,b] = min(abs(time/3600-23.1));
% DataLen = b-iD0; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% % [~,iD0] = min(abs(time/3600-22.1));
% % DataLen = 500; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% BeamList = [1 1 1 0];

% % % ARISTO 2017 RF04 - testing leg
% [~,iD0] = min(abs(time/3600-23));
% % [~,b] = min(abs(time/3600-23));
% DataLen = size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% % [~,iD0] = min(abs(time/3600-22.1));
% % DataLen = 500; %size(lams_spectra1,2)-iD0; %9241 size(lams_spectra1,2)
% BeamList = [1 1 1 0];

dt = mean(diff(time));

Tl2aTset = eye(3);  % LAMS to aircraft extra transformation matrix

timeWL = time(iD0:(iD0+DataLen-1));

ewx(isnan(ewx)) = 0;  % cases where water vapor pressure is nan confuse the algorithm.  Set them to zero.

% adjrateOffset = 2e-2; % 0.01
% adjrateCoef = 2e-5;
adjrate = 1e-2;  % learning rate of the coefficients

Ns = 1024;          % number of fft sample points
fs = 200e6;         % sample rate
dfs = fs/Ns;        % frequency resolution
lambda = 1560e-9;   % LAMS wavelength

% Slams = [(k1L).'; (k2L).'; (k3L).'; (k4L).'];

% % Half assed filter for bad INS data used on ARISTO2017 RF03
% iFill_INS = find(abs((roll-croll)-mean(roll-croll))>0.04);
% roll(iFill_INS) = roll(iFill_INS-1);
% pitch(iFill_INS) = pitch(iFill_INS-1);
% thdg(iFill_INS) = thdg(iFill_INS-1);

Tl2aVec = diffRPYvector(croll(:).',cpitch(:).',cthdg(:).',roll(:).',pitch(:).',thdg(:).');
Tl2a_Mean = nanmean(Tl2aVec(:,find(tas>50)),2);
Tl2a_Std = nanstd(Tl2aVec(:,find(tas>50)).').';
DeviationSNR = abs(Tl2aVec-Tl2a_Mean)./Tl2a_Std;
iUseMean_Tl2a = find(sqrt(sum(DeviationSNR.^2))>8 | isnan(sqrt(sum(DeviationSNR.^2))));
iUseData_Tl2a = find(sqrt(sum(DeviationSNR.^2))<8 & ~isnan(sqrt(sum(DeviationSNR.^2))));
% Tl2a_Mean_M = reshape(Tl2a_Mean,3,3);
% Tl2a_Mean_M = Tl2a_Mean_M/det(Tl2a_Mean_M);

Tl2aAvg = AvgTransformationMatrix(Tl2aVec(:,iUseData_Tl2a));
Tl2a_Mean_M = reshape(Tl2aAvg,3,3);


MaxTolerance = 12;  % absolute maximum deviation from TAS in m/s
BeamAngle = 30*pi/180;

memory_pts = 100;

Md = 28.9637;       % Mass of dry air in kg/kmol
Mw = 18.0153;       % Mass of water in kg/kmol
R0 = 8.314472e3;    % universal gas constant in J/kmol/K
Rd = R0/Md;         % gas constant for dry air
c_pd = 7*Rd/2;      % specific heat of dry air at constant pressure
c_vd = 5*Rd/2;      % specific heat of dry air at constant volume
c_pW = 4*R0;        % specific heat of water at constant pressure
c_vW = 3*R0;        % specific heat of water at constant volume
epsilon = Mw/Md;    % ratio of water weight to dry air
gamma_d = c_pd/c_vd;  % Ratio of dry specific heats

% tas_offset = 1.3;
% tas_std = 0.5;
% aa_offset = -0.02;%-0.02;
% aa_std = 0.01;
% ss_offset = 0.0024;%0;
% ss_std = 5e-3;

% Overall flight envelope definitions for aircraft flight parameters
ssenvel = 0.01;  % sideslip envelope
aaenvel = 0.02;  % angle of attack envelope
aamean =  0.03;  % mean angle of attack

% evelope definition for possible change in LOS velocity
da_std = 2*0.91;  % std of change in LOS speed between 1 sec data points


tas_std = 1.04;  % Estimated uncertainty in aircraft TAS estimate
aa_std = 0.0093; % Estimated uncertainty in aircraft angle of attack estimate
ss_std = 0.051; % Estimated uncertainty in aircraft sideslip estimate

% Difference between SND-500 and Honeywell
std_croll_d = std(croll(iD0-1+(1:DataLen))-roll(iD0-1+(1:DataLen)));
std_cpitch_d = std(cpitch(iD0-1+(1:DataLen))-pitch(iD0-1+(1:DataLen)));
std_cthdg_d = std(unwrap(cthdg(iD0-1+(1:DataLen)))-unwrap(thdg(iD0-1+(1:DataLen))));

% Uncertainty specifications for SDN-500 (from wind uncertainty document)
std_croll = 1e-3;
std_cpitch = 1e-3;
std_cthdg = 1.5e-3;

% Uncertainty specifications for Honeywell (from wind uncertainty document)
std_roll = 0.05*pi/180;
std_pitch = 0.05*pi/180;
std_thdg = 0.2*pi/180;

MinDataPts= 18;  % minimum error points needed to perform a correction to the coefficients matrix

if strcmp(Aircraft,'GV')
    % GV Initial Kalman coefficient defintions
    
    % Coefficients use ADIFR, BDIFR and QCF
    CoeffVector = zeros(7,1);
    CoeffVector(1) = 1.39; % TAS offset
    CoeffVector(2) = 0.99; % TAS multiplier
    CoeffVector(3) = 0.075; % angle of attack offset
    CoeffVector(4) = 0.304; % angle of attack multiplier
    CoeffVector(5) = 3.3642e-5;      % angle of attack machnumber-aoa multiplier
    CoeffVector(6) = 0.00093;    % side slip offset
    CoeffVector(7) = 0.20335; % side slip multiplier
    CoefLabel = {'TAS h_0','TAS h_1','\alpha f_0','\alpha f_1', '\alpha f_2','\beta g_0','\beta g_1'};

    % Coefficients use ADIFR, BDIFR and QCF
    % limits on coefficients
    xLim = zeros(size(CoeffVector,1),2);
    xLim(1,:) = [-20, 20];
    xLim(2,:) = [0.95, 1.15];
    xLim(3,:) = [0.051, 0.103];
    xLim(4,:) = [0.2511, 0.3945];
    xLim(5,:) = [-2e-4,2e-4];
    xLim(6,:) = [-0.0044,0.016];
    xLim(7,:) = [0.1477, 0.5908];  %[0.4, 1.1]
    aaMult = CoeffVector(4);
    ssMult = CoeffVector(7);
    tasMult = CoeffVector(2);

else
    % C-130 Initial Kalman coefficient defintions
    
    % ADIFR,BDIFR and QCF based coefficients
    CoeffVector = zeros(7,1);
    CoeffVector(1) = 2.7785; %-13.02;    % TAS offset
    CoeffVector(2) = 0.96; %1.0854;      % TAS multiplier
    CoeffVector(3) = 0.082656; %0.01;     % angle of attack offset
    CoeffVector(4) = 0.2117; %0.8731;     % angle of attack multiplier
    CoeffVector(5) = 0;      % angle of attack machnumber-aoa multiplier
    CoeffVector(6) = 0.019866;    % side slip offset
    CoeffVector(7) = 0.17928; %1.516; %0.652;     % side slip multiplier
    CoefLabel = {'TAS h_0','TAS h_1','\alpha f_0','\alpha f_1', '\alpha f_2','\beta g_0','\beta g_1'};

    % ADIFR,BDIFR and QCF based coefficients
    % limits on coefficients
    xLim = zeros(size(CoeffVector,1),2);
    xLim(1,:) = [-20, 20];
    xLim(2,:) = [0.95, 1.15];
    xLim(3,:) = [-0.011, 0.1725];
    xLim(4,:) = [0.092, 0.35];
    xLim(5,:) = [-2e-4,2e-4];
    xLim(6,:) = [8e-4,0.0532];
    xLim(7,:) = [0.0897, 0.3589];  %[0.0897, 0.3589]
end

% Build LAMS beam pointing matrix
BeamIndex = 1;
Slams = zeros(sum(BeamList),3);
if BeamList(1)
    Slams(BeamIndex,:) = k1L;
    BeamIndex = BeamIndex+1;
end
if BeamList(2)
    Slams(BeamIndex,:) = k2L;
    BeamIndex = BeamIndex+1;
end
if BeamList(3)
    Slams(BeamIndex,:) = k3L;
    BeamIndex = BeamIndex+1;
end
if BeamList(4)
    Slams(BeamIndex,:) = k4L;
end

% initialization of data vectors
vLAMSwl = zeros(3,DataLen);  % wavelet processing airspeed vector in LAMS coordinate frame
vLAMSac = zeros(3,DataLen);  % wavelet processing airspeed vector from LAMS in aircraft coordinate frame
vRD = zeros(3,DataLen);  % expected 3D velocities based on RD
vRDEst = zeros(3,DataLen);  % expected 3D velocities based on adaptive RD coeff
vLAMSglob = zeros(3,DataLen);  % wavelet processing airspeed vector from LAMS in global coordinate frame
vOmegaL = zeros(3,DataLen);  % correction due to aircraft rotation LAMS frame (uses SDN-500)
vOmegaL2 = zeros(3,DataLen);  % correction due to aircraft rotation LAMS frame (uses Honeywell)
vOmegaAC = zeros(3,DataLen);  % correction due to aircraft rotation Aircraft Frame (uses Honeywell)
vOmegaAC2 = zeros(3,DataLen);  % correction due to aircraft rotation Aircraft Frame (uses SDN-500)
vOmegaG = zeros(3,DataLen);  % correction due to aircraft rotation Global Frame
wLAMS = zeros(3,DataLen);  % 3D winds from LAMS
wEst = zeros(3,DataLen);  % 3D winds from Kalman Filtered RD
EstFlightParams = zeros(3,DataLen); % adaptively estimated flight parameters: [tas,ss,attack] from Kalman filter
EstErrors = zeros(3,DataLen); % difference between Kalman filtered radome and LAMS estimated flight parameters: [tas,ss,attack]
CoeffRecord = zeros(length(CoeffVector),DataLen);  % Kalman filter (sensitivity) coefficients
diffSS = zeros(3,DataLen);  % history of derivative (sensitivity) of LOS velocity to flight parameters
diffAA = zeros(3,DataLen);  % history of derivative (sensitivity) of LOS velocity to flight parameters
diffTAS = zeros(3,DataLen);  % history of derivative (sensitivity) of LOS velocity to flight parameters
dmach_ddp = zeros(1,DataLen);  % history of derivative (sensitivity) of mach number to pressure correction
dtas_mach = zeros(1,DataLen);  % history of derivative (sensitivity) of TAS to mach number

dp_p_corr_record = zeros(1,DataLen);  % history of pressure correction

% variance terms in the LAMS estimates
covLAMSwl= zeros(3,DataLen);                % variance of LAMS frame LAMS winds
covLAMSac = zeros(3,DataLen);               % variance of Aircraft frame LAMS winds
covLAMSLOS = zeros(sum(BeamList),DataLen);  % variance of LOS velocities from LAMS
covLAMSglob = zeros(sum(BeamList),DataLen); % variance of Global frame LAMS winds
var_tasLAMS = zeros(1,DataLen);             % variance of LAMS determined TAS
foundPk = zeros(sum(BeamList),DataLen);     % boolian indicating a peak was found in the beam spectrum. dims are # beams x length(time)
pPkrecord = zeros(sum(BeamList),DataLen);   % record of the selected peak probability

% Carries only the statistical portion of the covariance matrix
covLAMSwlStat= zeros(3,DataLen);                % variance of LAMS frame LAMS winds
covLAMSacStat = zeros(3,DataLen);               % variance of Aircraft frame LAMS winds
covLAMSLOSStat = zeros(sum(BeamList),DataLen);  % variance of LOS velocities from LAMS
covLAMSglobStat = zeros(sum(BeamList),DataLen); % variance of Global frame LAMS winds
covLAMSglobINS = zeros(sum(BeamList),DataLen); % variance of Global frame Including INS Error
covLAMSacINS = zeros(sum(BeamList),DataLen); % variance of Aircraft frame Including INS Error

errorKalman = zeros(sum(BeamList),DataLen);     % Kalman error in aircraft reference frame.  Includes both LAMS difference and statistical variation.

aDiff = zeros(sum(BeamList),DataLen);           % Difference between LAMS selected LOS velocities and Kalman predicted
aLAMSlist = zeros(sum(BeamList),DataLen);       % LAMS selected LOS velocities

NoOp = zeros(sum(BeamList),DataLen);

tasLAMS = zeros(1,DataLen);             % TAS determined from LAMS data
ssLAMS = zeros(1,DataLen);              % sideslip determined from LAMS data
attackLAMS = zeros(1,DataLen);          % angle of attack determined from LAMS data

tas_offset_array = zeros(1,DataLen);
ss_offset_array = zeros(1,DataLen);
aa_offset_array = zeros(1,DataLen);

tas_coef_array = zeros(1,DataLen);
ss_coef_array = zeros(1,DataLen);
aa_coef_array = zeros(1,DataLen);
tas_aa_coef_array = zeros(1,DataLen);

peakProbability = zeros(1,DataLen);
peakProbabilityAdj = zeros(1,DataLen);

% Kalman filtered radome parameters
tasEstv = zeros(1,DataLen);
ssEstv = zeros(1,DataLen);
aaEstv = zeros(1,DataLen);

aEstLast = zeros(sum(BeamList),1);  % initialize Kalman estimate projected onto the beams
estUncert = zeros(1,DataLen);       % uncertainty in Kalman estimate of projected LOS velocities

% lists of all possible peak locations
indexSelect = zeros(sum(BeamList),DataLen);  % index of selected peak for each beam (out of found peak list)
iSet = cell(sum(BeamList),DataLen);          % list of indicies into FFT spectrum of peak locations (out of 512 pt spectrum)
pSet = cell(sum(BeamList),DataLen);          % list of associated probabilities assosiated with each peak found in the FFT spectrum
spdSet = cell(sum(BeamList),DataLen);        % list of LOS speeds associated with each peak found in the FFT spectrum

mach_uc = zeros(1,DataLen);  % uncorrected mach number
mach_c = zeros(1,DataLen);  % corrected mach number
p_c = zeros(1,DataLen);  % corrected pressure

LR = adjrate;  % adjustable learning rate used in SVD Kalman filter with limits
LRcap = 0.1;  % maximum allowed learning rate in Kalman filter
LRrecord = zeros(1,DataLen);  % record of the SVD Kalman filter learning rate
xdim = length(CoeffVector);     % size of the Kalman state vector (number of coefficients)

WLspec1 = zeros(size(lams_spectra1,1),DataLen-1);   % Wavelet spectra for beam 1
WLspec2 = zeros(size(lams_spectra2,1),DataLen-1);   % Wavelet spectra for beam 2
WLspec3 = zeros(size(lams_spectra3,1),DataLen-1);   % Wavelet spectra for beam 3
% WLspec4 = zeros(size(lams_spectra4,1),DataLen-1); % Wavelet spectra for beam 4 

Lw1_spd = [];
Lw2_spd = [];
Lw3_spd = [];
Lw4_spd = [];

% Calculate correction to airspeeds due to lever arm length of LAMS from
% the aircraft center of mass.
% determine rotation vector of aircraft
Om = diff([unwrap(roll(:)) unwrap(pitch(:)) unwrap(thdg(:))])./diff(time(:));  % calculated from Honeywell data
Omc = diff([unwrap(croll(:)) unwrap(cpitch(:)) unwrap(cthdg(:))])./diff(time(:));  % Caculated from SDN-500 data

% Filter for large shifts in the data
Omf = Om;
Omf(abs(Omf) >= 0.1*pi/dt) = 0;

Omcf = Omc;
Omcf(abs(Omcf) >= 0.1*pi/dt) = 0;

radome_error = zeros(DataLen); % flag for nans in the radome data


% Iterate through point-by-point of flight data to find the optimal set of
% LAMS peaks in each time step.
for iprof = 1:DataLen

    DataIndex = iD0 +iprof-1;  % index into netcdf loaded data
    
    % Calculate transfomration matrices between coordinate frames
    Ta2g = RPYmat(roll(DataIndex),pitch(DataIndex),thdg(DataIndex)); % aircraft to global
    Tl2g = RPYmat(croll(DataIndex),cpitch(DataIndex),cthdg(DataIndex)); % LAMS to global
    % Check if INS seems to be working correctly.  If not use the average
    % transformation matrix between LAMS and Aircraft.
    if sum(DataIndex == iUseMean_Tl2a)
        Tl2a = Tl2aTset*Tl2a_Mean_M;
    else
        Tl2a = Tl2aTset*(Ta2g\Tl2g);  % LAMS to aircraft
    end
    
    % Compute velocity due to rotation of aircraft
    % where LAMS is on a lever arm relative to the aircraft center of mass
    vOmegaAC(:,iprof) = cross(Ta2g.'*(Omf(DataIndex,:).'),rL);  % Compute LAMS velocity error in aircraft frame using the Honeywell data
    vOmegaL(:,iprof) = cross(Tl2g.'*(Omcf(DataIndex,:).'),Tl2a.'*rL);  % Compute LAMS velocity error in LAMS frame using the SDN-500 data
    vOmegaG(:,iprof) = Tl2g*vOmegaL(:,iprof);
    if RotationCorrect
%         vRotK = -vOmegaAC(:,iprof); 
        vRotK = zeros(3,1);
    else
        vRotK = zeros(3,1);
    end
    
    vOmegaL2(:,iprof) = Tl2a.'*vOmegaAC(:,iprof); % Compute LAMS velocity error in LAMS frame using the Honeywell data 
    vOmegaAC2(:,iprof) = Tl2a*vOmegaL(:,iprof); % Compute LAMS velocity error in Aircraft frame using the SDN-500 data 
    
    % record the Kalman state vector coefficients used in this processing
    % step
    CoeffRecord(:,iprof) = CoeffVector;
    
    
    
    % Interpolate missing aircraft data (if it is missing)
    if isnan(qcf(DataIndex))
        % interpolate qcf
        radome_error(iprof) = 1;
        i1i = find(~isnan(qcf(1:DataIndex)),1,'last');
        i2i = find(~isnan(qcf(DataIndex:end)),1,'first')-1+DataIndex;
        qcf(DataIndex) = interp1([i1i,i2i],qcf([i1i,i2i]),DataIndex);
    end
    if isnan(bdifr(DataIndex))
        % Interpolate bdifr
        radome_error(iprof) = 1;
        i1i = find(~isnan(bdifr(1:DataIndex)),1,'last');
        i2i = find(~isnan(bdifr(DataIndex:end)),1,'first')-1+DataIndex;
        bdifr(DataIndex) = interp1([i1i,i2i],bdifr([i1i,i2i]),DataIndex);
    end
    if isnan(adifr(DataIndex))
        % Interpolate bdifr
        radome_error(iprof) = 1;
        i1i = find(~isnan(adifr(1:DataIndex)),1,'last');
        i2i = find(~isnan(adifr(DataIndex:end)),1,'first')-1+DataIndex;
        adifr(DataIndex) = interp1([i1i,i2i],adifr([i1i,i2i]),DataIndex);
    end
    
    % Estimate aircraft flight parameters from Kalman sensitivity
    % coefficients
    ssEst = [0,0,0,0,0,1,bdifr(DataIndex)/qcf(DataIndex)]*CoeffVector;  % Kalman adjusted sideslop
    aaEst = [0,0,1,adifr(DataIndex)/qcf(DataIndex),mach(DataIndex)*adifr(DataIndex)/qcf(DataIndex),0,0]*CoeffVector;  % Kalman adjusted angle of attack
    tasEst = [1,tas(DataIndex),0,0,0,0,0]*CoeffVector;  % Kalman adjusted TAS
   
    EstFlightParams(:,iprof) = [tasEst;ssEst;aaEst];  % Store the estimated flight parameters in an array
    
    % save Kalman filtered Radome Parameters
    tasEstv(iprof) = tasEst;
    aaEstv(iprof) = aaEst;
    ssEstv(iprof) = ssEst;

    % Compute flight derivatives (as a function of Kalman coefficients) for covariance matrix
    Sm = Slams*Tl2a.';  % Transform beam pointing projection to accept aircraft coordinate frame airflow input
    
    
    dSS = real(Sm*[0; sec(ssEst).^2; 0]*(tasEst));  % derivative of observed LOS speed with respect to side slip
    dAA = real(Sm*[0; 0; sec(aaEst).^2]*(tasEst));  % derivative of observed LOS speed with respect to angle of attack
    dTAS = real(Sm*[1; tan(ssEst);tan(aaEst)]);     % derivative of observed LOS speed with respect to TAS
    aEst = real(Sm*[1; tan(ssEst);tan(aaEst)]*(tasEst))+vRotK;  % Estimate LOS velocities for each beam based on Kalman filtered Radome data
    if iprof > 1
        daEst_dt = abs(aEst-aEstLast)*dfs*lambda/2/PkRepeatDerate;  % calculate time deriviative of Kalman estimated LOS velocities if this isn't the first data point
    else
        daEst_dt = zeros(sum(BeamList),1);  % if this is the first data point, set the LOS time derivative to zero
    end
    aEstLast = aEst;  % update the last LOS velocity measurement to be the current one
    
    % Compute the covariance matrix of the expected LOS observations based
    % on assumed uncertainty in sideslip, angle of attack and TAS
    CovA = real(dSS*ss_std.^2*dSS.'+dAA*aa_std.^2*dAA.'+dTAS*tas_std.^2*dTAS.');
    
    % Save the derivatives of the flight parameters
    diffSS(:,iprof) = dSS;
    diffAA(:,iprof) = dAA;
    diffTAS(:,iprof) = dTAS;
    
    % Calculate covariance matrix for an overall aircraft flight envelope
    CovAenvel = dSS*ssenvel.^2*dSS.'+dAA*aaenvel.^2*dAA.'+dTAS*tas_std.^2*dTAS.';
    aEnv = Sm*[1; tan(0);tan(aamean)]*(tasEst);  % set the mean value of the flight envelope
    

    %%% Find spectral peaks in the recorded LAMS spectra
    spectHist = max([2,DataIndex-30]);
    if BeamList(1)
        lams_spectra1(1,DataIndex) = lams_spectra1(2,DataIndex);  % remove zero frequency term and replace with adjacent bin
        [b1pks,b1spd,ib1,pib1,NoOp(1,iprof),WLspec] = PeakFindIndex_wavelet(aWL,lams_spectra1(:,DataIndex),BeamAngle,tas(DataIndex),MaxTolerance,2);  % find the peaks in the spectrum using wavelet analysis
        WLspec1(:,iprof) = WLspec;  % store the wavelet spectra for aWL (wavelet width in fft points) = 2
        
        % determine if any of the found peaks need to be derated for
        % appearing repeatedly in the same location (possible noise spurs)
        Lw1_new = ones(length(ib1{1}),1);
        for ipk = 1:length(ib1{1})
           if sum(round(ib1{1}(ipk)) == Lw1_spd)
               Lw1_new(ipk) = prod(Lw1(round(ib1{1}(ipk)) == Lw1_spd))*nanmax([nanmin([PkRepeatDerate/daEst_dt(1),1]),PkRepeatDerate]);
           end
        end
        Lw1_spd = [0; round(ib1{1})];  % extra term to account for estimation "peak" in case a suitible observed peak is not found
        Lw1 = [NoPkProb; Lw1_new];  % extra term to account for estimation "peak"
    end
    if BeamList(2)
        lams_spectra2(1,DataIndex) = lams_spectra2(2,DataIndex);%0.5*lams_spectra2(1,DataIndex);
        [b2pks,b2spd,ib2,pib2,NoOp(2,iprof),WLspec] = PeakFindIndex_wavelet(aWL,lams_spectra2(:,DataIndex),0,tas(DataIndex),MaxTolerance,2);
        WLspec2(:,iprof) = WLspec;
        
        % determine if any of the found peaks need to be derated for
        % appearing repeatedly in the same location (possible noise spurs)
        Lw2_new = ones(length(ib2{1}),1);
        for ipk = 1:length(ib2{1})
           if sum(round(ib2{1}(ipk)) == Lw2_spd)
               Lw2_new(ipk) = prod(Lw2(round(ib2{1}(ipk)) == Lw2_spd))*nanmax([nanmin([PkRepeatDerate/daEst_dt(2),1]),PkRepeatDerate]);
           end
        end
        Lw2_spd = [0; round(ib2{1})];  % extra term to account for estimation "peak"
        Lw2 = [NoPkProb; Lw2_new];  % extra term to account for estimation "peak"
    end
    if BeamList(3)
        lams_spectra3(1,DataIndex) = lams_spectra3(2,DataIndex);%0.5*lams_spectra2(1,DataIndex);
%         spec_input = lams_spectra3(:,DataIndex) - mean(lams_spectra3(:,spectHist:DataIndex-1),2);
%         [b3pks,b3spd,ib3,pib3,NoOp(3,iprof),WLspec] = PeakFindIndex_wavelet(aWL,spec_input,BeamAngle,tas(DataIndex),MaxTolerance,2);
        [b3pks,b3spd,ib3,pib3,NoOp(3,iprof),WLspec] = PeakFindIndex_wavelet(aWL,lams_spectra3(:,DataIndex),BeamAngle,tas(DataIndex),MaxTolerance,2);
        WLspec3(:,iprof) = WLspec;
        
        % determine if any of the found peaks need to be derated for
        % appearing repeatedly in the same location (possible noise spurs)
        Lw3_new = ones(length(ib3{1}),1);
        for ipk = 1:length(ib3{1})
           if sum(round(ib3{1}(ipk)) == Lw3_spd)
               Lw3_new(ipk) = prod(Lw3(round(ib3{1}(ipk)) == Lw3_spd))*nanmax([nanmin([PkRepeatDerate/daEst_dt(3),1]),PkRepeatDerate]);
           end
        end
        Lw3_spd = [0; round(ib3{1})];  % extra term to account for estimation "peak"
        Lw3 = [NoPkProb; Lw3_new];  % extra term to account for estimation "peak"
    end
    if BeamList(4)
        lams_spectra4(1,DataIndex) = lams_spectra4(2,DataIndex);%0.5*lams_spectra2(1,DataIndex);
%         spec_input = lams_spectra4(:,DataIndex) - mean(lams_spectra4(:,spectHist:DataIndex-1),2);
%         [b4pks,b4spd,ib4,pib4,NoOp(4,iprof),WLspec] = PeakFindIndex_wavelet(aWL,spec_input,BeamAngle,tas(DataIndex),MaxTolerance,2);
        [b4pks,b4spd,ib4,pib4,NoOp(4,iprof),WLspec] = PeakFindIndex_wavelet(aWL,lams_spectra4(:,DataIndex),BeamAngle,tas(DataIndex),MaxTolerance,2);
        WLspec4(:,iprof) = WLspec;
        
        Lw4_new = ones(length(ib4{1}),1);
        for ipk = 1:length(ib4{1})
           if sum(round(ib4{1}(ipk)) == Lw4_spd)
               Lw4_new(ipk) = prod(Lw4(round(ib4{1}(ipk)) == Lw4_spd))*nanmax([nanmin([PkRepeatDerate/daEst_dt(4),1]),PkRepeatDerate]);
           end
        end
        Lw4_spd = [0; round(ib4{1})];  % extra term to account for estimation "peak"
        Lw4 = [NoPkProb; Lw4_new];  % extra term to account for estimation "peak"
    end
    %%% end peak finding rountine

    foundPk(:,iprof) = 1;
    % update beam uncertainties based on peak selections
    % if no peaks were found in a particular beam, use the Kalman
    % determined peak and assign the associated high uncertainty to that
    % observation.
    PkLen = 1;
    if BeamList(1)
        if isempty(b1spd{1})
            b1spd{1} = aEst(1);
            ib1{1} = nan;   
            std1 = stdNoPk;  % store error associated with no peak found (usually 10m/s)
            foundPk(1,iprof) = 0;
        else
            b1spd{1} = [aEst(1); b1spd{1}];
            ib1{1} = [nan; ib1{1}];
            std1 = max(std(b1spd{1}),dfs*lambda/2);
        end
        PkLen = PkLen*length(b1spd{1});
    end
    if BeamList(2)
        if isempty(b2spd{1})
            b2spd{1} = aEst(2);
            ib2{1} = nan;
            std2 = stdNoPk;
            foundPk(2,iprof) = 0;
        else
            b2spd{1} = [aEst(2); b2spd{1}];
            ib2{1} = [nan; ib2{1}];
            std2 = max(std(b2spd{1}),dfs*lambda/2);
        end
        PkLen = PkLen*length(b2spd{1});
    end
    if BeamList(3)
        if isempty(b3spd{1})
            b3spd{1} = aEst(3);
            ib3{1} = nan;
            std3 = stdNoPk;
            foundPk(3,iprof) = 0;
        else
            b3spd{1} = [aEst(3); b3spd{1}];
            ib3{1} = [nan; ib3{1}];
            std3 = max(std(b3spd{1}),dfs*lambda/2);
        end
        PkLen = PkLen*length(b3spd{1});
    end
    if BeamList(4)
        if isempty(b4spd{1})
            b4spd{1} = aEst(4);
            ib4{1} = nan;
            std4 = stdNoPk;
            foundPk(4,iprof) = 0;
        else
            b4spd{1} = [aEst(4); b4spd{1}];
            ib4{1} = [nan; ib4{1}];
            std4 = max(std(b4spd{1}),dfs*lambda/2);
        end
        PkLen = PkLen*length(b4spd{1});
    end
    
    % normCoef = sqrt(det(CovA)/pi^4);
    % Normalization coefficient for Normal PDFs.  Set to 1 because
    % otherewise we are dealing with floating point numbers at the edge of
    % their precision.
    normCoef = 1;
    
    % Evaluate every peak combination found in the wavelet spectral
    % analysis to determine what combination is most likely.
    ProbPks = zeros(sum(BeamList)+1,PkLen);
    indSet = zeros(sum(BeamList),PkLen);
    iProbPks = 1;
    for i1 = 1:length(b1spd{1})
        for i2 = 1:length(b2spd{1})
            for i3 = 1:length(b3spd{1})
%                 for i4 = 1:length(b4spd{1})
                    % Three beam version
                    xSpd = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3)]-aEst;  % Evaluate difference between this peak combination and expected (Kalman) velocities
                    xEnv = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3)]-aEnv;  % limits based on aircraft flight parameters
%                     % Four beam version
%                     xSpd = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4)]-aEst;
%                     xEnv = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4)]-aEnv;  % limits based on aircraft flight parameters

%                     normCoef = ((i1==1)*NoPkProb+(i1~=1)*1.0)*((i2==1)*NoPkProb+(i2~=1)*1.0)*((i3==1)*NoPkProb+(i3~=1)*1.0);
                    normCoef = Lw1(i1)*Lw2(i2)*Lw3(i3);  % evaluate penalties associated with this peak combination (e.g. repeating peak or Kalman determined peak)
%                      normCoef = Lw1(i1)*Lw2(i2)*Lw3(i3)*Lw4(i4);
                    
                    if iprof > 1
                        % time derivative evaluations (requires history of
                        % data to calculate)
                        
%                         dxSpd = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4)]-aLAMS-dadtEst;
%                         ProbPks(:,iProbPks) = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4);normCoef*(exp(-xSpd.'*CovA*xSpd)*exp(-dxSpd.'*CovdAdt*dxSpd))]; 
%                         dxSpd = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4)]-(aLAMS+mean(diff([aLAMSlist(:,max([1,iprof-3]):iprof-1)-aDiff(:,max([1,iprof-3]):iprof-1) aEst].')).');
                        
                        % Change in LOS speed of each beam for a given
                        % point selection
                        dxSpd = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3)]-(aLAMS+mean(diff([aLAMSlist(:,max([1,iprof-3]):iprof-1)-aDiff(:,max([1,iprof-3]):iprof-1) aEst].')/dt).');
                        
                        % Check if the covariance matrix of the Kalman
                        % solution is singular (or close to it).
                        if det(CovA+probAdj+SolnCov) < 1e-10
                            % Singular covariance matrix, use
                            % psuedo-inverse to project lower dimensional
                            % Gaussian into beam dimension space (3 or 4)
                            CovAInv = pinv(CovA+probAdj+SolnCov);
                        else
                            % Non-signular covariance.  Compute Gaussian
                            % pdf in the full beam dimension space.
                            CovAInv = inv(CovA+probAdj+SolnCov);
                        end
                        % Calculate the probability that this particular
                        % peak combination is the right one by assuming all
                        % error terms can be described using Normal PDFs.
                        ProbPks(:,iProbPks) = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);normCoef*exp(-xSpd.'*CovAInv*xSpd).*exp(-xEnv.'*inv(CovAenvel)*xEnv).*exp(-dxSpd.'*((1/(da_std*max([1,estUncert(iprof-1)/50])))).^2*dxSpd)];
                    else
                        % first data point, so no time derivative
                        % information can be used.
                        dxSpd = zeros(4,1);
                        if det(CovA) < 1e-10
                            CovAInv = pinv(CovA);
                        else
                            CovAInv = inv(CovA);
                        end
%                         ProbPks(:,iProbPks) = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);b4spd{1}(i4);normCoef*exp(-xSpd.'*(CovA)*xSpd).*exp(-xEnv.'*(CovAenvel)*xEnv)];
                        ProbPks(:,iProbPks) = [b1spd{1}(i1);b2spd{1}(i2);b3spd{1}(i3);normCoef*exp(-xSpd.'*CovAInv*xSpd).*exp(-xEnv.'*inv(CovAenvel)*xEnv)]; 
                    end
                    % Store the indices to the peak in each beam for later
                    % use.
%                     indSet(:,iProbPks) = [i1;i2;i3;i4];  % four beam
                    indSet(:,iProbPks) = [i1;i2;i3];  % three beam
                    iProbPks = iProbPks+1;  % increment index for the next peak combination
%                 end
            end
        end
    end
    
    %find the most likely peak combination
    [maxProbPks,iaLAMS] = max(ProbPks(sum(BeamList)+1,:));
    
    % calculate probability weighted LOS beam velocity covariance matrix
    pNorm = ProbPks(sum(BeamList)+1,:)/sum(ProbPks(sum(BeamList)+1,:));
    SolnCov = zeros(sum(BeamList));
    for iC1 = 1:sum(BeamList)
        for iC2 = 1:sum(BeamList)
            SolnCov(iC1,iC2) = sum(pNorm.*ProbPks(iC1,:).*ProbPks(iC2,:))-sum(pNorm.*ProbPks(iC1,:))*sum(pNorm.*ProbPks(iC2,:));
        end
    end
    
    aLAMS = ProbPks(1:sum(BeamList),iaLAMS);  % LOS speeds determined from LAMS peak selection
    aDiff(:,iprof) = aLAMS-aEst;    % difference between Kalman filtered radome estimated and LAMS estimated LOS speeds
    aLAMSlist(:,iprof) = aLAMS;     % store LAMS LOS speeds in an array
    
    peakProbability(iprof) = maxProbPks;  % store the evaluated probability of the selected peak combination
    
    % Compute covariance matrix scaling factor based on success of
    % predicting the previous 6 points in the flight
    if iprof > 6
        probScale = exp(-(0:5)).*prod(foundPk(:,iprof:-1:iprof-5),1)./(1+da_std*(0:5));
        if sum(probScale) ~= 0
            probScale = probScale/sum(probScale);
            % use an exponential average, but also reduce the probabilities
            % based on the amount of time since the measurement
%             probAdj = sum(probScale.*peakProbability(iprof:-1:iprof-5)); 

              probAdj = sqrt(sum(probScale.*(aDiff(:,iprof:-1:iprof-5).^2),2));
%               probAdj = diag(probAdj.^2);
              probAdj = probAdj*probAdj.';
        else
%             probAdj = 0.1;
            if BeamList(4)
                probAdj = diag([std1, std2, std3, std4]);
            else
                probAdj = diag([std1, std2, std3]);
            end
        end
    else
        probAdj = 0.1;
    end
    
%     peakProbabilityAdj(iprof) = probAdj;
    
    % update beam uncertainties based on peak selections
    aCovDiag = [];
    if BeamList(1)
        if ~isempty(pib1{1}) && indSet(1,iaLAMS) ~= 1
            std1 = max(std(b1spd{1}),dfs*lambda/2*ResFactor)/exp(-(1-pib1{1}(indSet(1,iaLAMS)-1))^2);
            pPkrecord(1,iprof) = pib1{1}(indSet(1,iaLAMS)-1);
        else
            foundPk(1,iprof) = 0;
            std1 = stdNoPk;
            SolnCov(1,1) = SolnCov(1,1)+std1^2;
        end
        
        SolnCov(1,1) = SolnCov(1,1) + (dfs*lambda/2*ResFactor)^2;
        
        indexSelect(1,iprof) = ib1{1}(indSet(1,iaLAMS));
        iSet{1,iprof} = ib1{1};
        pSet{1,iprof} = pib1{1};
        spdSet{1,iprof} = b1spd{1};
        std1 = sqrt(std1^2 + (aLAMS(1)-aEst(1)).^2/(CovA(1,1)));
        
        aCovDiag = [aCovDiag 1./std1^2];
    end
    if BeamList(2)
        if ~isempty(pib2{1}) && indSet(2,iaLAMS) ~= 1
            std2 = max(std(b2spd{1}),dfs*lambda/2*ResFactor)/exp(-(1-pib2{1}(indSet(2,iaLAMS)-1))^2);
            pPkrecord(2,iprof) = pib2{1}(indSet(2,iaLAMS)-1);
        else
            foundPk(2,iprof) = 0;
            std2 = stdNoPk;
            SolnCov(2,2) = SolnCov(2,2)+std2^2;
        end
        
        SolnCov(2,2) = SolnCov(2,2) + (dfs*lambda/2*ResFactor)^2;
        
        indexSelect(2,iprof) = ib2{1}(indSet(2,iaLAMS));
        iSet{2,iprof} = ib2{1};
        pSet{2,iprof} = pib2{1};
        spdSet{2,iprof} = b2spd{1};
        std2 = sqrt(std2^2 + (aLAMS(2)-aEst(2)).^2/(CovA(2,2)));
        
        aCovDiag = [aCovDiag 1./std2^2];
    end
    if BeamList(3)
        if ~isempty(pib3{1}) && indSet(3,iaLAMS) ~= 1
            std3 = max(std(b3spd{1}),dfs*lambda/2*ResFactor)/exp(-(1-pib3{1}(indSet(3,iaLAMS)-1))^2);
            pPkrecord(3,iprof) = pib3{1}(indSet(3,iaLAMS)-1);
        else
            foundPk(3,iprof) = 0;
            std3 = stdNoPk;
            SolnCov(3,3) = SolnCov(3,3)+std3^2;
        end
        
        SolnCov(3,3) = SolnCov(3,3) + (dfs*lambda/2*ResFactor)^2;
        
        indexSelect(3,iprof) = ib3{1}(indSet(3,iaLAMS));
        iSet{3,iprof} = ib3{1};
        pSet{3,iprof} = pib3{1};
        spdSet{3,iprof} = b3spd{1};
        std3 = sqrt(std3^2 + (aLAMS(3)-aEst(3)).^2/(CovA(3,3)));
        
        aCovDiag = [aCovDiag 1./std3^2];
    end
    if BeamList(4)
        if ~isempty(pib4{1}) && indSet(4,iaLAMS) ~= 1
            std4 = max(std(b4spd{1}),dfs*lambda/2*ResFactor)/exp(-(1-pib4{1}(indSet(4,iaLAMS)-1))^2);
            pPkrecord(4,iprof) = pib4{1}(indSet(4,iaLAMS)-1);
        else
            foundPk(4,iprof) = 0;
            std4 = stdNoPk;
            SolnCov(4,4) = SolnCov(4,4)+std4^2;
        end
        
        SolnCov(4,4) = SolnCov(4,4) + (dfs*lambda/2*ResFactor)^2;
        
        indexSelect(4,iprof) = ib4{1}(indSet(4,iaLAMS));
        iSet{4,iprof} = ib4{1};
        pSet{4,iprof} = pib4{1};
        spdSet{4,iprof} = b4spd{1};
        std4 = sqrt(std4^2 + (aLAMS(4)-aEst(4)).^2/(CovA(4,4)));
        
        aCovDiag = [aCovDiag 1./std4^2];
    end
    
    
%     aCovInv = diag(aCovDiag);
    aCovInv = inv(SolnCov+probAdj);
    
    if iprof > 1
        estUncert(iprof) = sum(ProbPks(sum(BeamList)+1,:))/ProbPks(sum(BeamList)+1,iaLAMS)+estUncert(iprof-1)*0.9;
    else
        estUncert(iprof) = sum(ProbPks(sum(BeamList)+1,:))/ProbPks(sum(BeamList)+1,iaLAMS);
    end
    
    % calculate airflow parameters and wind parameters
    
    % 
    % assume TAS is only x component of airflow measurement
    vRD(:,iprof) = [1; tan(ssrd(DataIndex));tan(attackrd(DataIndex))]*(tas(DataIndex));  
%     % assume TAS is the magnitude of the total airflow measurement
%     vRD(:,iprof) = [1; tan(ssrd(DataIndex));tan(attackrd(DataIndex))]*(tas(DataIndex))/sqrt(1+tan(ssrd(DataIndex)).^2+tan(attackrd(DataIndex)).^2);

    vRDEst(:,iprof) = [1; tan(ssEst);tan(aaEst)]*(tasEst);

    vInvMat = ((((Slams.')*aCovInv)*Slams)\((Slams.')*aCovInv));  % velocity inversion matrix.  For overdefined measurements (4 beams) it performs an uncertainty weighted estimate.
    vLAMSwl(:,iprof) =  vInvMat*aLAMS; % velocities in LAMS frame
%     vLAMSwl(:,iprof) = Slams\aLAMS;
    vLAMSac(:,iprof) = Tl2a*vLAMSwl(:,iprof);  % velocities in aircraft frame
    
    vLAMSglob(:,iprof) = Tl2g*(vLAMSwl(:,iprof)-vOmegaL(:,iprof));  % LAMS velocities in global frame
    wLAMS(:,iprof) = vLAMSglob(:,iprof) -[cvns(DataIndex); cvew(DataIndex); cvspd(DataIndex)];  % 3D winds from LAMS
    wEst(:,iprof) = Ta2g*vRDEst(:,iprof) - [vnsc(DataIndex); vewc(DataIndex); -gvspd(DataIndex)];  % 3D winds from Kalman filtered radome
    
%     % Error Terms
%     cov_vLAMSwl = vInvMat*(SolnCov+probAdj)*vInvMat.';
%     var_vLAMSwl(:,iprof) = diag(cov_vLAMSwl);
%     var_vLAMSac(:,iprof) = diag(Tl2a*vInvMat*(SolnCov+probAdj)*(Tl2a*vInvMat).');
%     var_vLAMSglob(:,iprof) = diag(Tl2g*vInvMat*(SolnCov+probAdj)*(Tl2g*vInvMat).');
    
    % calculate flight parameters from LAMS (use TAS = x component of
    % aircraft velocity)
%     tasLAMS(iprof) = vLAMSac(1,iprof);
%     ssLAMS(iprof) = atan(vLAMSac(2,iprof)./tasLAMS(iprof));
%     attackLAMS(iprof) = atan(vLAMSac(3,iprof)./tasLAMS(iprof));

    % calculate flight parameters from LAMS (use TAS = |v|)
    tasLAMS(iprof) = sqrt(sum(vLAMSac(:,iprof).^2));
    ssLAMS(iprof) = asin(vLAMSac(2,iprof)./tasLAMS(iprof));
    attackLAMS(iprof) = asin(vLAMSac(3,iprof)./tasLAMS(iprof));
    
    

%     covLAMSwlMat = vInvMat*diag(1./aCovDiag)*vInvMat.';
    covLAMSwlMat = vInvMat*(SolnCov+probAdj)*vInvMat.';  % covariance matrix of LAMS air motion estimate
    covLAMSwlMatStat = vInvMat*(SolnCov)*vInvMat.';  % covariance matrix of LAMS air motion estimate without Kalman error (statistical only error)
    
    covLAMSwl(:,iprof) = diag(covLAMSwlMat);  % store full LAMS air motion covariance matrix
    covLAMSac(:,iprof) = diag(Tl2a*covLAMSwlMat*Tl2a.');  % store diagonals of aircraft frame LAMS air motion covariance matrix
    covLAMSLOS(:,iprof) = diag(SolnCov+probAdj);           % store diagonals of LAMS LOS covariance matrix
    covLAMSglob(:,iprof) = diag(Tl2g*vInvMat*(SolnCov+probAdj)*(Tl2g*vInvMat).');  % store diagonals of global frame LAMS covariance matrix
    
    covLAMSwlStat(:,iprof) = diag(covLAMSwlMat);
    covLAMSacStat(:,iprof) = diag(Tl2a*covLAMSwlMatStat*Tl2a.');
    covLAMSLOSStat(:,iprof) = diag(SolnCov);
    covLAMSglobStat(:,iprof) = diag(Tl2g*vInvMat*(SolnCov)*(Tl2g*vInvMat).');
%     covLAMSLOS(:,iprof) = 1./aCovDiag;

    [dTdroll,dTdpitch,dTdyaw] = diffRPYmat(croll(DataIndex),cpitch(DataIndex),cthdg(DataIndex));
    covLAMSglobINS(:,iprof) = covLAMSglobStat(:,iprof) + diag( ...
        dTdroll*vLAMSwl(:,iprof)*std_croll^2*(dTdroll*vLAMSwl(:,iprof)).' + ...
        dTdpitch*vLAMSwl(:,iprof)*std_cpitch^2*(dTdpitch*vLAMSwl(:,iprof)).' + ...
        dTdyaw*vLAMSwl(:,iprof)*std_cthdg^2*(dTdyaw*vLAMSwl(:,iprof)).');
    
    [dTdroll_ac,dTdpitch_ac,dTdyaw_ac] = diffRPYmat(croll(DataIndex)-roll(DataIndex),cpitch(DataIndex)-pitch(DataIndex),cthdg(DataIndex)-thdg(DataIndex));
    covLAMSacINS(:,iprof) = covLAMSacStat(:,iprof) + diag( ...
        dTdroll_ac*vLAMSwl(:,iprof)*(std_croll^2+std_roll^2)*(dTdroll_ac*vLAMSwl(:,iprof)).' + ...
        dTdpitch_ac*vLAMSwl(:,iprof)*(std_cpitch^2+std_pitch^2)*(dTdpitch_ac*vLAMSwl(:,iprof)).' + ...
        dTdyaw_ac*vLAMSwl(:,iprof)*(std_cthdg^2+std_thdg^2)*(dTdyaw_ac*vLAMSwl(:,iprof)).');
    
    
    % estimated error in the LAMS measurement based on the error in Kalman
    % measurement and propagation of statistical error (treated separately)
    errorKalman(:,iprof) = sqrt(covLAMSacStat(:,iprof))+0.5*abs(vLAMSac(:,iprof)-vRDEst(:,iprof));
    
    dtas_dvLAMSwl = vLAMSwl(:,iprof)/tasLAMS(iprof);
%     var_tasLAMS(iprof) = dtas_dvLAMSwl.'*covLAMSwlMat*dtas_dvLAMSwl;
    var_tasLAMS(iprof) = covLAMSwlStat(1,iprof);
    
    % update calibration parameters
    iMem1 = max(iD0,DataIndex-memory_pts);
    iMem2 = max(1,iprof-memory_pts);
    
%     % adapt the offset parameters based on measurements
    if sum(foundPk(:,iprof)) && tas(DataIndex) > 50 && radome_error(iprof)==0
        EstErrors(:,iprof) = -[tasEst-tasLAMS(iprof); ssEst-ssLAMS(iprof); aaEst-attackLAMS(iprof)];
        
        % find the point where we have enough data for a defined or
        % overdefined error analysis
        
        % Allows peaks to be in only one channel
%         PkHist = sum(foundPk(:,1:iprof));
%         PkHist = cumsum(PkHist,'reverse');
%         iHist = find(PkHist >= MinDataPts,1,'last');
        
        % Requires even distribution of peaks in channels to provide a
        % correction
        PkHist = (radome_error(1:iprof)==0).*foundPk(:,1:iprof);
        PkHist = cumsum(PkHist,2,'reverse');
%         iHist = find(prod(PkHist >= MinDataPts/sum(BeamList),1),1,'last');  % Require points be distributed between all 3 beams
        iHist = find(sum(PkHist,1)>= MinDataPts,1,'last');  % Allow points to be distributed in any way in the 3 beams

        if ~isempty(iHist) && (iprof-iHist) < 60
            DataiHist = iHist + iD0-1;  % index into the start of the flight record used to adjust the sensitivity coefficients

            rowSel = foundPk(:,iHist:iprof);
            rowSel = rowSel(:);

%             WeightedError = aDiff(:,iHist:iprof);     % Use for pinv
            weights = foundPk(:,iHist:iprof)./covLAMSLOSStat(:,iHist:iprof);  %sqrt(covLAMSLOSStat
            WeightedError = real(aDiff(:,iHist:iprof).*weights);  % Use for SVD
            WeightedError = WeightedError(:);
            weights = weights(:);
            Jacobian = zeros(sum(BeamList)*(iprof-iHist+1),length(CoeffVector));  % Coefficient list: tas_offset, tascoeff, tas_mach_coeff,tas_aa_coeff, aa_offset,aacoeff, aa_mach_coeff, ss_offset,sscoeff

            % Linear offsets and coefficeints
            Jacobian(:,1) = reshape(diffTAS(:,iHist:iprof)*1,[],1);
            Jacobian(:,2) = reshape(diffTAS(:,iHist:iprof).*tas(DataiHist:DataIndex).',[],1);
            Jacobian(:,3) = reshape(diffAA(:,iHist:iprof)*1,[],1);
            Jacobian(:,4) = reshape(diffAA(:,iHist:iprof).*(adifr(DataiHist:DataIndex)./qcf(DataiHist:DataIndex)).',[],1);
            Jacobian(:,5) = reshape(diffAA(:,iHist:iprof).*(adifr(DataiHist:DataIndex)./qcf(DataiHist:DataIndex).*mach(DataiHist:DataIndex)).',[],1);
            Jacobian(:,6) = reshape(diffSS(:,iHist:iprof).*1,[],1);
            Jacobian(:,7) = reshape(diffSS(:,iHist:iprof).*(bdifr(DataiHist:DataIndex)./qcf(DataiHist:DataIndex)).',[],1);

            % remove components where a peak wasn't found
            Jacobian(rowSel==0,:) = [];
            WeightedError(rowSel==0) = [];
            weights(rowSel==0) = [];
            CovError = covLAMSLOS(:,iHist:iprof);
            CovError = CovError(:);
            CovError(rowSel==0) = [];
            CovError = diag(1./CovError);
            

            
            % Use SVD to eliminate ineffective search directions in the Jacobian
            xmask = ones(size(CoeffVector));
            xm2 = ones(size(Jacobian));
            JFx = Jacobian.*weights;
            Fx = -WeightedError;
            [uJ,sJ,vJ]=svd(xm2.*JFx);
            if size(sJ,2) == 1
                sJv = sJ(1,1);
            else
                sJv = diag(sJ);
            end
            LR = min([LR*1.1,LRcap]);
            mJv=(Fx'*uJ(:,1:min(xdim,size(uJ,2))))./sJv';
            mJv(abs(sJv) < 1e-4*max(max(abs(sJv)))) = 0;
            mJ = diag(mJv);
            vJ = vJ(:,1:size(mJ,1));  % Rescale v if system is underdefined
            CoeffNew = CoeffVector - LR*sum(vJ*mJ,2).*xmask;

            % check if the solutions are out side of their allowable
            % limits.
            % if they are, remove the dominant principle components driving
            % us over the limit and reduce the learning rate.  Do this
            % until the solution is within the parameter space limits.
            LR0 = LR;
            adjIter = 0;
            lltest = 1;
            ultest = 1;
            while lltest || ultest
                while lltest
                    xll = xLim(:,1)>CoeffNew;
                    if sum(xll) > 0
                        if sum(diag(mJ)~=0) > 1
                            [~,bll] = max((vJ(xll,:)*mJ).');
                            mJ(bll,bll) = 0;
                            LR = LR0*0.969*exp(-adjIter^2/(11.5)^2);
                            CoeffNew = CoeffVector - LR*sum(vJ*mJ,2).*xmask;
                            ultest = 1;
                        else
                            LR = LR0*0.969*exp(-adjIter^2/(11.5)^2);
                            CoeffNew = CoeffVector - LR*sum(vJ*mJ,2).*xmask;
                        end
                        adjIter = adjIter +1;
                    else
                        lltest = 0;
                    end
                end


                while ultest
                    xul = xLim(:,2)<CoeffNew;
                    if sum(xul) > 0
                        if sum(diag(mJ)~=0) > 1
                            [~,bul] = min((vJ(xul,:)*mJ).');
                            mJ(bul,bul) = 0;
                            LR = LR0*0.969*exp(-adjIter^2/(11.5)^2);
                            CoeffNew = CoeffVector - LR*sum(vJ*mJ,2).*xmask;
                            lltest = 1;
                        else
                            LR = LR0*0.969*exp(-adjIter^2/(11.5)^2);
                            CoeffNew = CoeffVector - LR*sum(vJ*mJ,2).*xmask;
                        end
                        adjIter = adjIter +1;
                    else
                        ultest = 0;
                    end
                end
            end
            CoeffVector = real(CoeffNew);
            LRrecord(iprof) = LR;
        end
        

    else
        EstErrors(:,iprof) = nan*ones(3,1);
    end

end

% apply correction for lever arm rotation to observations
if RotationCorrect
   vLAMSac = vLAMSac-vOmegaAC; 
   vLAMSwl = vLAMSwl-vOmegaL;
%    vLAMSglob = vLAMSglob-vOmegaG; %Tl2g*vOmegaL(:,iprof);  % already performed in the processing loop
   vRDEst = vRDEst-vOmegaAC2;  
end

foundNan = foundPk;
foundNan(foundPk==0) = nan;

% Difference in total airspeed is an indicator for the type of error in the
% LAMS/radome comparision.  If the error is a transformation, then there
% should be no difference in total airspeed between the two measurements.
% The difference only appears when the actual measurements have error in
% them.
errorTotalAirspeed = (sqrt(sum(vRD.^2,1))-sqrt(sum(vLAMSwl.^2,1)));
errorTotalAirspeed2 = vRD(1,:)-sqrt(sum(vLAMSwl.^2,1));
errorTotalAirspeed3 = sqrt(sum(vRD.^2,1))-aLAMSlist(2,:);
errorTotalAirspeed4 = vRD(1,:)-aLAMSlist(2,:);

