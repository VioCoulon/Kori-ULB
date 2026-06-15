function KoriMeltStandalone

global_path='Z:\Vio\COUPLING\';
% global_path='Z:\Vio\OceanIce\FORCING_DATA\';

global_path = 'Z:\Vio\COUPLING\';
path_data   = [global_path,'Kori/input_Kori/'];
OCN_clim    = 'NN'; % ISMIP6 or NEMO or NN

resolution = 8;
init_name  = ['Bedmachine',int2str(resolution),'km']; % ALREADY CROPPED!!!

ctr.meltfunc    = 23; % melt scheme -- 3: PICO - 23: QUAD mean Ant slope (Burgard22) - 24: QUAD local slope (Burgard22)
gammaT          = 3e-5; % To adapt according to chosen melt scheme
mixedgamma      = 0; 
path_mixedgamma = [path_data,'kmixed_Antslope_extra.nc'];%[path_data,'kmixed_locslope_extra.nc'];
ctr.meltfac     = 1; 
ctr.gammaTplume = 0; % needs to be defined if meltfunc=5

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

ctr.shelf=1;
ctr.SSA=2;
ctr.calving=2;

ctr.dt=1;
ctr.nsteps=1;
ctr.diagnostic=1;

% OCEAN CLIM
if isequal(OCN_clim,'ISMIP6')
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

load([path_data,'ZBextended_',int2str(resolution),'km'],'ZB'); % make sure we use the right input basins - here Zwally
save(init_name,'ZB','-append')

addpath('C:\Users\Vio\OneDrive - Université Libre de Bruxelles\KoriDev\') % If not on cluster, make sure to use KoriCoupling Kori version
addpath('C:\Users\Vio\OneDrive - Université Libre de Bruxelles\KoriDev\subroutines')

KoriModel(init_name,['MELT_',int2str(ctr.meltfunc),'_',int2str(resolution),'km'],ctr,fc);

end