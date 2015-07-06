clear all; close all; clc;

% Load in LDAS point output in netcdf format
% Save as a .mat file of needed variables for each calender year
% Then use Combine_matlab_yearly_files.m to combine all calender year matlab files together

%% Site info %% IMPORTANT CHECK HERE!! %
timeoffset= -8; 
% Entire time period (will download one year at a time)
STIME_all     = datenum(1979,01,01,0,0,0);
ETIME_all     = datenum(2015,05,27,0,0,0);
TIME_OUT_all  = time_builder(STIME_all,ETIME_all,1);

% Paths
filesin  = '/usr/lusers/nicway/civil/NLDAS/netcdffiles';
filesout = '/usr/lusers/nicway/civil/NLDAS/matlabfiles/';

% Split period into calender years
Cyears = 1981:2015;

% Grid cell to extract
lat1  = 47.4380;
lon1 = -121.4380;
% Elevation from grid cell (very very close to NLDAS actual center)
%  29 180 47.4375 -121.4375   1126.2670
NLDAS_elev = 1126.2670; % meters
% http://ldas.gsfc.nasa.gov/nldas/NLDASelevation.php


%% Do work

for cY = Cyears 

% Get date array for current year
I_time_c = find(TIME_OUT_all(:,1)==cY);
TIME_OUT = TIME_OUT_all(I_time_c,:);
Nout     = size(TIME_OUT,1);


cd(filesin)
list = dir(['NLDAS_FORA0125_H.A' num2str(cY) '*']);
Nfiles = length(list);


LDAS_Data = nan(Nout,7);
VariableNLDAS     = {'SW down','LW down','Temp','Specific humidity','Precip','Wind Speed','Pressure'}; 
UnitsNLDAS        = {'W/m^2','W/m^2','K','kg/kg','kg/m^2','m/s','Pa'};

tic

for cf = 1:Nfiles;
    
    ncfile = list(cf).name;
    
    % Get Time
    ctime =  datenum(nc_attget(ncfile, 'TMP_110_HTGL', 'initial_time'));
    
    % Check if we want this time
    Icf = find(TIME_OUT(:,7) == ctime);
    
    if ~isempty(Icf)
        
        % Get lat and long vectors
        Lats = nc_varget(ncfile,'lat_110');
        Lons = nc_varget(ncfile,'lon_110');
        
        % Find the grid cell we want
        tol = 0.00001; 
        Ilat = find( abs(Lats - lat1) < tol );
        Ilon = find( abs(Lons - lon1) < tol );
        
        if isempty(Ilat) | isempty(Ilon)
            error('lat or lon not found')
        end
        
        gridcell = [Ilat,Ilon] - 1; % Subtract 1 to go from 1 index to 0 based index

        
        % Average wind speed
        V = nc_varget(ncfile,'V_GRD_110_HTGL',gridcell,[1 1]); % V
        U = nc_varget(ncfile,'U_GRD_110_HTGL',gridcell,[1 1]); % U
        W = sqrt(V^2 + U^2);

        % Get Data from cell we want
        LDAS_Data(Icf,1) = nc_varget(ncfile,'TMP_110_HTGL',gridcell,[1 1]);          % Temp K
        LDAS_Data(Icf,2) = nc_varget(ncfile,'SPF_H_110_HTGL',gridcell,[1 1]);         % Specific humidiy kg/kg at 2m
        LDAS_Data(Icf,3) = nc_varget(ncfile,'A_PCP_110_SFC_acc1h',gridcell,[1 1]);   % Precip kg/m^2
        LDAS_Data(Icf,4) = W;                                                          % Wind Speed m/s
        LDAS_Data(Icf,5) = nc_varget(ncfile,'DSWRF_110_SFC',gridcell,[1 1]);        % SW down W/m^2
        LDAS_Data(Icf,6) = nc_varget(ncfile,'DLWRF_110_SFC',gridcell,[1 1]);        % LW down W/m^2
        LDAS_Data(Icf,7) = nc_varget(ncfile,'PRES_110_SFC',gridcell,[1 1]);         % Pressure Pa
    end
    
    clear ncfile ctime Icf V U W
end
toc

cd(filesout)

%% Convert from GMT to PST (-8)

sprintf('Saving year %i\n',cY)
TIME_OUT = time_shift(TIME_OUT,timeoffset);
fileoutCY = ['NLDAS_data_CY_' num2str(cY)];
eval(['save ' fileoutCY ' TIME_OUT LDAS_Data VariableNLDAS UnitsNLDAS lat1 lon1 NLDAS_elev'])

clear TIME_OUT LDAS_Data VariableNLDAS UnitsNLDAS

end

disp('Finished')

return

