function var = NCloaddata(varName,ncfilename,SampleTime)

ncid = netcdf.open(ncfilename, 'NC_NOWRITE');

if ~isnan(sum(SampleTime))
    varid_time = netcdf.inqVarID(ncid,'Time');
    time = double(netcdf.getVar(ncid,varid_time));
end

varid = netcdf.inqVarID(ncid,varName);
var = double(netcdf.getVar(ncid,varid));
netcdf.close(ncid);

% convert angles in degrees to radians
try
    no_units = 0;
    varUnits = ncreadatt(ncfilename,varName,'units');
    if strcmp(varUnits, 'degree_T') || strcmp(varUnits, 'degree') && ~strcmp(varName,'ATX')
        var = var*pi/180;
        var = unwrap(var(:));
    end
catch
    no_units = 1;
end

if ~isnan(sum(SampleTime))
    varSR = numel(var)/numel(time);

    if varSR ~= 1
        timeVar = time(1)+(1:(length(time)*varSR))*mean(diff(time))/varSR;
    else
        timeVar = time;
    end

    if numel(var) ~= length(SampleTime) && numel(SampleTime) > 1
        var = interp1(timeVar,var(:),SampleTime);
    end
end
    



var(var == -32767) = nan;
var = var(:);
if no_units == 0
    if strcmp(varUnits, 'degree_T') || strcmp(varUnits, 'degree')
        var = mod(var+pi,2*pi)-pi;
    %     var = atan(sin(var)./cos(var));
    end
end