

ncfilenameB = 'ARISTO2017rf03.nc';  % LAMS Data file
ncfilenameO = 'ARISTO2017rf03.nc';  % Aircraft Data File
UHSAS_Suffix = '';  % Suffix for UHSAS variable names if it is in the instrument payload
addpath('/scr/raf_data/ARISTO2017');
Aircraft = 'GV';  % Use if this is a GV flight

% ncfilenameB = 'ARISTO2016rf06.nc';  % LAMS Data file
% ncfilenameO = 'ARISTO2016rf06.nc';  % Aircraft Data File
% UHSAS_Suffix = '_LWI';  % Suffix for UHSAS variable names if it is in the instrument payload
% addpath('/scr/raf_data/ARISTO2016/bin1');
% Aircraft = 'C-130'; % Use if this is a C-130 flight

SaveFigs = 0;  % automatically save figures
WriteNetCDF = 0;
SaveDirectory = '/scr/sci/mhayman/LAMS/';

FileTag = '_winds';

LAMS_LoadFlight

% set process indices (full flight)
% iD0 = 1;% 2000
% [~,iD0] = min(abs(time/3600-19.065));  % ARISTO2016rf06
[~,iD0] = min(abs(time/3600-20.08));  % ARISTO2017 rf04
DataLen = size(lams_spectra1,2)-iD0;
BeamList = [1 1 1 0];  

LAMS_ProcessFlight;

LAMS_PlotResults;

WriteWaveletNetcdf;