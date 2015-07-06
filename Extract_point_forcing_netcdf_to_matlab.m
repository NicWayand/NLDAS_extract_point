clear all; close all; clc;

% Load in LDAS point output in netcdf format
% Save as a .mat file of needed variables for each water year

% LDAS download instructions

% cd('G:\ncar\d1\data\LDAS\raw\SNQ_2014_Oct_May')
%cd('G:\ncar\d1\data\LDAS\raw\SNQ_continuous')
cd('/usr/lusers/nicway/civil/NLDAS/netcdffiles')

%% Site info %% IMPORTANT CHECK HERE!! %
timeoffset= -8; 
% STIME     = datenum(2012,09,30);
% ETIME     = datenum(2013,06,30);
% STIME     = datenum(2013,06,30);
% ETIME     = datenum(2013,11,2);
% % STIME     = datenum(2013,10,1,8,0,0);
% % ETIME     = datenum(2014,5,17,12,0,0);
%STIME     = datenum(2012,09,30,0,0,0);
STIME     = datenum(1979,01,01,0,0,0);
ETIME     = datenum(2015,05,27,0,0,0);
TIME_OUT  = time_builder(STIME,ETIME,1);
Nout      = size(TIME_OUT,1);

% % % Old cell (too far to east!!!! Do not use!!!)
% % % Grid cell to extract
% % lat1  = 47.4380;
% % lon1 = -121.3130;
% % % Elevation from grid cell (very very close to NLDAS actual center)
% % %  30 180 47.4375 -121.3125   1346.9290 % meters
% % % http://ldas.gsfc.nasa.gov/nldas/NLDASelevation.php

% Closest to SNQ - Nic - 6/11/2015
% Grid cell to extract
lat1  = 47.4380;
lon1 = -121.4380;
% Elevation from grid cell (very very close to NLDAS actual center)
%  29 180 47.4375 -121.4375   1126.2670
NLDAS_elev = 1126.2670; % meters
SNQ_elev = 921; % meters
% http://ldas.gsfc.nasa.gov/nldas/NLDASelevation.php

list = dir('NLDAS_FORA0125_H*');
Nfiles = length(list);

LDAS_Data = nan(Nout,7);
VariableNLDAS     = {'SW down','LW down','Temp','Specific humidity','Precip','Wind Speed','Pressure'}; 
UnitsNLDAS        = {'W/m^2','W/m^2','K','kg/kg','kg/m^2','m/s','Pa'};

tic
% Failed on 22098
% save mat_temp_dump.mat

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

        % Get Data from cell we want (center of 6 grid cells over SNQ)
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

cd('/usr/lusers/nicway/civil/NLDAS/matlabfiles/')

%% Convert from GMT to PST (-8)

TIME_OUT = time_shift(TIME_OUT,timeoffset);

% Uncomment below and change name to save
disp('Uncomment below and change name to save')
% save NLDAS_SNQ_ALL_raw.mat TIME_OUT LDAS_Data VariableNLDAS UnitsNLDAS lat1 lon1 NLDAS_elev SNQ_elev

return





% % % % return
% % % % clear all; close all; clc;
% % % % 
% % % % load('G:\ncar\d1\data\LDAS\raw\SNQ_continuous_out\NLDAS_SNQ_ALL_raw.mat')
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %% NLDAS values are instaneous values, except for Precipitation which is the sum of previous hour!
% % % % 
% % % % % Disaggregate SW and LW down from hourly to 30 min (to match forcing)
% % % % [TIME_30min T_30min]      = disaggFLUX(TIME_OUT, LDAS_Data(:,1)); %#ok<*ASGLU>
% % % % [TIME_30min SH_30min]     = disaggFLUX(TIME_OUT, LDAS_Data(:,2));
% % % % Precip_30min              = disaggMET(TIME_OUT(:,7),LDAS_Data(:,3),get_dt(TIME_OUT),get_dt(TIME_30min),TIME_30min(:,7),1, 1,0);
% % % % [TIME_30min W_30min]      = disaggFLUX(TIME_OUT, LDAS_Data(:,4));
% % % % [TIME_30min SW_30min]     = disaggFLUX(TIME_OUT, LDAS_Data(:,5));
% % % % SW2                       = SHIN_INST(TIME_OUT,LDAS_Data(:,5),TIME_30min,47.4249,-121.4138,8,'END');
% % % % [TIME_30min LW_30min]     = disaggFLUX(TIME_OUT, LDAS_Data(:,6));
% % % % [TIME_30min Press_30min]  = disaggFLUX(TIME_OUT, LDAS_Data(:,7));
% % % % 
% % % % % Convert stuff
% % % % Press_30min = Pa_to_hPa(Press_30min); 
% % % % UnitsNLDAS{7} = 'hPa';
% % % % 
% % % % %% Adjust/correct NLDAS values to SNQ elevation
% % % % 
% % % % %% Shift Temperature from MAU elevation to Station gridcell avg elevation
% % % % % T_Lapse = -5.5; % Annual from Minder. Calculate if we can see a seasonal correction here.
% % % % % % Significant with three years?
% % % % % SNQ_elev = 921; % meters
% % % % % T_shift       = T_Lapse/1000*(SNQ_elev-NLDAS_elev);
% % % % % fprintf('Tshift is %f C from %f m to %f m\n',T_shift,NLDAS_elev,SNQ_elev)
% % % % % T_30min = T_30min + T_shift;
% % % % 
% % % % % % Shift Longwave (based on Marks(?) study that found -29 W/m^2 per km)
% % % % % LW_lapse = -29; % W/m^2 per km)
% % % % % LW_shift           = LW_lapse/1000*(Pix_elev-NLDAS_elev);
% % % % % fprintf('LW shift is %f W/m^2 from %f m to %f m\n',LW_shift,NLDAS_elev,Pix_elev)
% % % % % LW_MAU = LW_MAU + LW_shift;
% % % % 
% % % % 
% % % % 
% % % % %% apply corrections calculated with Calculate_NLDAS_biases_to_observation_point.m
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % % % Plot to check
% % % % % figure(1)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,5),'-k*')
% % % % % plot(TIME_30min(:,7),SW_30min,'-b*')
% % % % % plot(TIME_30min(:,7),SW2,'-r*')
% % % % % legend('NLDAS 1hr INST','orig disagg','New disagg')
% % % % % tlabel
% % % % % 
% % % % % [EL, AZ, SOLDIST, HA, DEC] = SolarGeometry_v2(TIME_30min,47.4249,-121.4138,8);
% % % % % 
% % % % % 
% % % % % figure(2)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,6),'-k*')
% % % % % plot(TIME_30min(:,7),LW_30min,'-b*')
% % % % % tlabel
% % % % % 
% % % % % figure(3)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,1),'-k*')
% % % % % plot(TIME_30min(:,7),T_30min,'-b*')
% % % % % tlabel
% % % % % 
% % % % % figure(4)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,2),'-k*')
% % % % % plot(TIME_30min(:,7),SH_30min,'-b*')
% % % % % tlabel
% % % % % 
% % % % % figure(5)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,3),'-k*')
% % % % % plot(TIME_30min(:,7),Precip_30min,'-b*')
% % % % % tlabel
% % % % % 
% % % % % if nansum(Precip_30min) ~= nansum(LDAS_Data(:,3))
% % % % %     error('Precip disagg changed total')
% % % % % end
% % % % % 
% % % % % figure(6)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,4),'-k*')
% % % % % plot(TIME_30min(:,7),W_30min,'-b*')
% % % % % tlabel
% % % % % 
% % % % % figure(7)
% % % % % hold on
% % % % % plot(TIME_OUT(:,7),LDAS_Data(:,7),'-k*')
% % % % % plot(TIME_30min(:,7),Press_30min,'-b*')
% % % % % tlabel
% % % % 
% % % % 
% % % % % Format for input using Import_X_2_Master_Format.m
% % % % data_NLDAS = [SW2 LW_30min T_30min SH_30min Precip_30min W_30min Press_30min];
% % % % 
% % % % % Save
% % % % save 'G:\ncar\d1\data\LDAS\matlabF\SNQ_WY_2013_2015.mat' data_NLDAS TIME_30min VariableNLDAS UnitsNLDAS
% % % % 
% % % % return





