function KoriMeltStandalone

    
clear;
close all;


% Paths.
addpath /home/daniel/models/Kori-ULB/subroutines/;
addpath /home/daniel/models/Kori-ULB/;


global_path = '/home/daniel/models/Kori-ULB/ice_data/ismip7/ismip7-antarctic-ocean-forcing';
%path_param  = [global_path,'/parameterisations/ocean/ocean_modelling_data/'];


% Mathiot_NEMO_cold_m.nc      melt_warm_target_term3.nc         Naughten_FESOM_MMM_cold_m.nc   Timmermann_FESOM_cold_v2_S.nc
%Mathiot_NEMO_cold_v2_S.nc   Naughten_FESOM_ACCESS_cold_m.nc   Naughten_FESOM_MMM_cold_S.nc   Timmermann_FESOM_cold_v2_TF.nc
%Mathiot_NEMO_cold_v2_TF.nc  Naughten_FESOM_ACCESS_cold_S.nc   Naughten_FESOM_MMM_cold_TF.nc  Timmermann_FESOM_cold_v2_T.nc
%Mathiot_NEMO_cold_v2_T.nc   Naughten_FESOM_ACCESS_cold_TF.nc  Naughten_FESOM_MMM_cold_T.nc   Timmermann_FESOM_warm_m.nc
%Mathiot_NEMO_warm_m.nc      Naughten_FESOM_ACCESS_cold_T.nc   Naughten_FESOM_MMM_warm_m.nc   Timmermann_FESOM_warm_v2_S.nc
%Mathiot_NEMO_warm_v2_S.nc   Naughten_FESOM_ACCESS_warm_m.nc   Naughten_FESOM_MMM_warm_S.nc   Timmermann_FESOM_warm_v2_TF.nc
%Mathiot_NEMO_warm_v2_TF.nc  Naughten_FESOM_ACCESS_warm_S.nc   Naughten_FESOM_MMM_warm_TF.nc  Timmermann_FESOM_warm_v2_T.nc
%Mathiot_NEMO_warm_v2_T.nc   Naughten_FESOM_ACCESS_warm_TF.nc  Naughten_FESOM_MMM_warm_T.nc
%melt_cold_target_term3.nc   Naughten_FESOM_ACCESS_warm_T.nc   Timmermann_FESOM_cold_m.nc

OCN_clim = 'ISMIP7'; % ISMIP6 or NEMO or NN


resolution = 8;
init_name  = ['Bedmachine',int2str(resolution),'km_v3_RACMO11km_Stal2021']; % ALREADY CROPPED!!!

ctr.runmode     = 3;
ctr.meltfunc    = 3; % melt scheme -- 3: PICO - 23: QUAD mean Ant slope (Burgard22) - 24: QUAD local slope (Burgard22)
ctr.C           = 1e6;
gammaT          = 3e-5; % To adapt according to chosen melt scheme
mixedgamma      = 0; 
%path_mixedgamma = [path_data,'kmixed_Antslope_extra.nc'];%[path_data,'kmixed_locslope_extra.nc'];
ctr.meltfac     = 1; 
ctr.gammaTplume = 0; % needs to be defined if meltfunc=5

% Control.
ctr.shelf      = 1;
ctr.SSA        = 2;
ctr.calving    = 2;
ctr.dt         = 1;
ctr.nsteps     = 1;
ctr.diagnostic = 1;



if mixedgamma==0
    ctr.gammaT=gammaT;
else
    gamma_field=ncread(path_mixedgamma,'coeff_k'); % updated field computed by Clara (B24)
    ctr.gammaT=gamma_field';
end

if resolution==16
    ctr.imax=351;
    ctr.jmax=351;
    ctr.delta=16.e3;
elseif resolution==8
    ctr.imax=701;
    ctr.jmax=701;
    ctr.delta=8.e3;
end


% Define range of parameters.
if ctr.meltfunc == 3
    gammaT_0 = 1e-5;
    gammaT_f = 1e-4;
    n_gammaT = 2;

    C_0 = 1e5;
    C_f = 2e6;
    n_C = 2;

    values_1 = linspace(C_0, C_f, n_C);
    values_2 = linspace(gammaT_0, gammaT_f, n_gammaT);

end



% OCEAN CLIM
if isequal(OCN_clim,'ISMIP7')

    path_to  = [global_path,'/obs/ocean/climatology/zhou_annual_06_nov/thetao/v4/'];
    path_so  = [global_path,'/obs/ocean/climatology/zhou_annual_06_nov/so/v4/'];

    file_to = 'thetao_AIS_obs_ocean_climatology_zhou_annual_06_nov_v4_1972-2024.nc';
    file_so = 'so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v4_1972-2024.nc';

    full_to = [path_to, file_to];
    full_so = [path_so, file_so];

    % Dimensions definition.
    info = ncinfo(full_so,'so');
    {info.Dimensions.Name}

    to = ncread(full_to, 'thetao');
    so = ncread(full_so, 'so');
    z  = ncread(full_so, 'z');

    size(to)

    % Crop domain.
    nISMIP = 761;
    nf = 701;
    dn = 0.5 * (n0 - nf);

    To = to(dn+1:end-dn, dn+1:end-dn, :);
    So = so(dn+1:end-dn, dn+1:end-dn, :);

    save(init_name,'To','So','-append')

    % Depths of the 30 layers, it should be provided with the climatology.
    % This is what tells Kori what depth is every layer for the 
    % interpolation at the ice shelf draft.
    fc.z = z;

elseif isequal(OCN_clim,'ISMIP6')
    load ([path_data,'ISMIP6_OCEAN_OBS_CLIMATOLOGY_1995-2017_',int2str(resolution),'km.mat'],'z')
    fc.z=z;
elseif isequal(OCN_clim,'NEMO')
    z=ncread([path_data,'clim-nemo_1982-2013.nc'],'deptht');
    z=double(z);
    fc.z=-z;
elseif isequal(OCN_clim,'NN')
    z=ncread([path_data,'clim-nemo_1982-2013.nc'],'deptht');
    z=double(z);
    fc.z=-z;
end

% Load data for ISMIP6 melt param
if ctr.meltfunc==9
    load([path_data,'ocean_param_',int2str(resolution),'km.mat']);
    fc.deltaT_basin=deltaT_basin;
    fc.basinNumber=basinNumber;
elseif ctr.meltfunc==91
    load([path_data,'ocean_param_slope_',int2str(resolution),'km.mat']);
    fc.deltaT_basin=deltaT_basin;
    fc.basinNumber=basinNumber;
end


if isequal(OCN_clim,'ISMIP6')
    load([path_data,'ISMIP6_OCEAN_OBS_CLIMATOLOGY_1995-2017_',int2str(resolution),'km.mat'],'To','So')
    save(init_name,'To','So','-append')
elseif isequal(OCN_clim,'NEMO')
    To=ncread([path_data,'clim-nemo_1982-2013.nc'],'thetao');
    To=permute(To,[2 1 3]);
    So=ncread([path_data,'clim-nemo_1982-2013.nc'],'so');
    So=permute(So,[2 1 3]);
    save(init_name,'To','So','-append')
elseif isequal(OCN_clim,'NN')
    To=ncread([path_data,'clim-nn_1982-2013.nc'],'thetao');
    To=permute(To,[2 1 3]);
    So=ncread([path_data,'clim-nn_1982-2013.nc'],'so');
    So=permute(So,[2 1 3]);
    save(init_name,'To','So','-append')
end

%load([path_data,'ZBextended_',int2str(resolution),'km'],'ZB'); % make sure we use the right input basins - here Zwally
%save(init_name,'ZB','-append')



path_out = '/home/daniel/models/Kori-ULB/output/calibration/pico/';
mkdir(path_out); 

%KoriModel(init_name,['MELT_',int2str(ctr.meltfunc),'_',int2str(resolution),'km'],ctr,fc);

file = [path_out,'MELT_',int2str(ctr.meltfunc),'_',int2str(resolution),'km','_'];

% Loop over range of parameters for a given parametrization choice.
for i = 1:length(values_1)
    for j = 1:length(values_2)

        ctr.C      = values_1(i);
        ctr.gammaT = values_2(j);

        val_1 = 1e-6*values_1(i);
        val_2 = 1e5*values_2(j);

        if val_1 < 10
            num_1 = sprintf('0%.0f', val_1);
        else
            num_1 = sprintf('%.0f', val_1);
        end

        if val_2 < 10
            num_2 = sprintf('0%.0f', val_2);
        else
            num_2 = sprintf('%.0f', val_2);
        end

        % Name of each permutation.
        name = ['C', num_1, '_', 'gamma', num_2];

        % Run Kori diagnostic exp.
        KoriModel(init_name, [file,name], ctr, fc);
        
        % Load Melt from source file
        S = load([file,name,'_toto'], 'Melt');

        % Store under dynamic field name
        Melt = zeros(nISMIP, nISMIP);
        melt_all.(name) = S.Melt;

    end
end

% Save ensemble of melt rates together.
%file_all = 'MELT_',int2str(ctr.meltfunc),'_',int2str(resolution),'km','_','melt_all.mat';
%file = [path_out,'MELT_',int2str(ctr.meltfunc),'_',int2str(resolution),'km','_'];
save([file,'melt_all.mat'], 'melt_all');

end