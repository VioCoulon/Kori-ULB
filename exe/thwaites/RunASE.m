function RunASE

% Script to run ASE at 2 km resolution.
% Kori-ULB v0.9

clear;
close all;

addpath /home/daniel/models/Kori-ULB/subroutines/;
addpath /home/daniel/models/Kori-ULB/;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First SIA inversion (As)
ctr.imax     = 472;
ctr.jmax     = 400;
ctr.delta    = 2.e3;
ctr.nsteps   = 4001; % 4001
ctr.dt       = 5;
ctr.Ao       = 5.0e-17; % 5.0e-17
ctr.m        = 3;
ctr.Asin     = zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.calving  = 2;
ctr.Tinit    = 1;    % Initialization of temperature field from semi-analytical steady-state temperature solution
ctr.Tcalc    = 2;    % Calculate temperature field and thermomechanical coupling, , i.e. A = f (T)
ctr.basin    = 1;    % Run the model for a specific basin.
ctr.inverse  = 1;  % optimization of basal sliding coefficients As for the grounded ice sheet with fixed grounding line position.
ctr.runmode  = 3;  % 1: graphics; 3: no graphics
%KoriModelAll('ASE2km','INITA_NON',ctr);
%KoriModel('ASE2km','INITA_NON',ctr);

path = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/ice_data/eta1e7/ground_melt_1/';

name_1 = 'ASE2km';
name_2 = 'INIT_SIA';
input  = strcat(path, name_1);
output = strcat(path, name_2);

%KoriModel(input, output, ctr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second SSA inversion (As, MeltInv) - short run
ctr.inverse      = 2;        % Optimization of basal sliding coefficients As for the grounded ice sheet and sub-shelf melt/accretion for floating ice shelves.
ctr.meltfunc     = 1;        % Beckmann and Goosse (2003) with linear dependency on the thermal forcing, following (de Boert al., 2015). 
ctr.GroundedMelt = 1;        % necessary for basins!! Vio: 0. Frank: 1.
ctr.shelf        = 1;        % Ice shelves are considered.
ctr.SSA          = 2;
ctr.nsteps       = 101;
ctr.dt           = 0.01;
ctr.Tinit        = 0;        % SET 1 FOR INITIALIZATION ! Initial temperature field read from input file or when not available kept constant at values of surface temperature
ctr.Tinv         = 10;       % Time interval between updates in the optimization scheme
ctr.TinvMelt     = 0.01;
ctr.HinvMelt     = 10;


name_1 = 'INIT_SIA';
name_2 = 'INIT_SSA_1';
input  = strcat(path, name_1);
output = strcat(path, name_2);

%KoriModel(input, output, ctr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second SSA inversion (As, MeltInv) - short run
ctr.nsteps   = 40001; % 10000, 40000
ctr.dt       = 0.1; % 0.1
ctr.TinvMelt = 5;
ctr.Tinv     = 5;
ctr.Hinv     = 200;

name_1 = 'INIT_SSA_1';   % INITB_NON
name_2 = 'INIT_SSA_2';  % INIT_NON_1
input  = strcat(path, name_1);
output = strcat(path, name_2);

%KoriModel(input, output, ctr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat longer run to stabilze
ctr.nsteps = 40001;  % Frank: 40000
%ctr.dt     = 0.1;   % 0.1

name_1 = 'INIT_SSA_2';  % INIT_NON_1
name_2 = 'INIT_SSA_3';  % INIT_NON_2
input  = strcat(path, name_1);
output = strcat(path, name_2);

%KoriModel(input, output, ctr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing run
ctr.inverse  = 0;
%ctr.calving = 2;          % ice shelves do not extend further than initially
ctr.SSA      = 2;
ctr.m        = 5.0;
ctr.calving  = 4;          % 5, to apply LSF function.
ctr.dt       = 0.05;        % 0.02, 0.1, 0.05, 0.025
ctr.nsteps   = 20001;       % 5001, 20001
ctr.meltfunc = 3;          % PICO
ctr.gammaT   = 2.0e-5;     % 2.0e4, 1.0e-4, 0.25e-4
ctr.meltfac  = 5;          % Factor multiplying sub-shelf melt.


ctr.timeslice = 1;
ctr.snapshot  = 50;

ctr.tforcing = 20.0

% Test stochastic forcing.
%ctr.stochastic = 0;

% Test damage with dynamic grain size.
ctr.bassis_reg = 0;

% Graphics.
ctr.runmode  = 1; 

% Define paths.
%path_1      = [parent_path, 'Initialization/eta1e7/ground_melt_1/'];
%path_2      = [parent_path, 'Deterministic/eta1e7/ground_melt_1/calv2/deter_gamma1e-2/'];
parent_path = '/home/daniel/models/Kori-ULB/output/Thwaites/';
%path_1      = [parent_path, 'Initialization/eta1e7/ground_melt_0/']; % 'Initialization/eta1e7/ground_melt_0/'
path_1      = [parent_path, 'Initialization/eta1e7/ground_melt_0/basal_friction/coulomb/m_05/']; % 'Initialization/eta1e7/ground_melt_0/'
path_2      = [parent_path, 'Deterministic/test/'];


init_name = 'INIT_SSA_3';             % INIT_SSA_3
out_name  = 'deter_bassis_reg1';        % deter_gamma1e-2
path_in   = strcat(path_1, init_name);
path_out  = strcat(path_2, out_name);

KoriModel(path_in, path_out, ctr);

% Make 2 runs: with and without stochastic forcing.
% Chech if we should limit the values of melt.
%KoriModel('INIT_NON_END','ASEfor_m1',ctr);

%KoriModel('INIT_NON_END_2','ASEfor_m1',ctr);

% Ensure the ocean temperatures for the basin are read.
%load('ASE2km','To');
%save('INIT_NON_3','To','-append');

%KoriModel('INIT_NON_3','ASEfor_m1_deter_gamma5e-4',ctr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ctr.damage=1;
% KoriModel('INITB_NON','ASEfor_m2',ctr);


end


