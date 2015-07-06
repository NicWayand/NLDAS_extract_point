clear all; close all; clc;

% Combine multiple matlab files of NLDAS yearly data

%% Site info %% IMPORTANT CHECK HERE!! %
STIME_all     = datenum(1979,01,01,0,0,0);
ETIME_all     = datenum(2015,05,27,0,0,0);
TIME_OUT_all  = time_builder(STIME_all,ETIME_all,1);
Nout          = size(TIME_OUT_all,1);
TIME_OUT_all = datenum_round_off(TIME_OUT_all,'minute');
LDAS_Data_all = nan(Nout,7);

% Paths
filesin = '/usr/lusers/nicway/civil/NLDAS/matlabfiles/';
fileoutALL = 'SNQ_1979_2015.mat';

%% Do work
cd(filesin)
list = dir('*NLDAS*.mat');


for cY = 1:length(list) 

    mtfile = list(cY).name;
    load(mtfile)
    
    TIME_OUT = datenum_round_off(TIME_OUT,'minute');
    
    [C,A,B] = intersect(TIME_OUT_all(:,7),TIME_OUT(:,7));
    
    LDAS_Data_all(A,:) = LDAS_Data(B,:);
    
    clear TIME_OUT LDAS_DATA
    
end

TIME_OUT = TIME_OUT_all;
NLDAS_Data = LDAS_Data_all;

% Save
eval(['save ' fileoutALL ' TIME_OUT NLDAS_Data VariableNLDAS UnitsNLDAS lat1 lon1 NLDAS_elev SNQ_elev'])


disp('Finished')

return

