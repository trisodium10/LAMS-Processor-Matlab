function var = NCloadhist(varName,ncfilename)
% var = NCloadhist(varName,ncfilename)
% load 2D histogram data from netCDF
% varName - string with variable name
% ncfilename - string for netCDF filename

ncid = netcdf.open(ncfilename, 'NC_NOWRITE');

varid = netcdf.inqVarID(ncid,varName);
var = double(netcdf.getVar(ncid,varid));
var = reshape(var,[size(var,1) size(var,2)*size(var,3)]);

netcdf.close(ncid);