function RunThwaites

% ASE initialisation at 2km
% Kori-ULB

clear;
close all;

%% First SIA inversion (As)
ctr.imax=472;
ctr.jmax=400;
ctr.delta=2.e3;
ctr.nsteps=8001;
ctr.dt=5;
ctr.Ao=5.0e-17; % 5.0e-17
ctr.m=3;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.calving=2; % fix calving front to current position
ctr.Tinit=1;
ctr.Tcalc=2;
ctr.basin=1;
ctr.inverse=1;
ctr.runmode=1;

if ctr.runmode==1
    addpath('C:\Users\Vio\OneDrive - Université Libre de Bruxelles\KoriDev\') % If not on cluster, make sure to use latest Kori version
    addpath('C:\Users\Vio\OneDrive - Université Libre de Bruxelles\KoriDev\subroutines')
end

%KoriModel('ASE2km','ASE2km_initsia',ctr);

%% Second SSA inversion (As, MeltInv) - short run

ctr.inverse=2;
ctr.meltfunc=1;
ctr.GroundedMelt=0; 
ctr.shelf=1;
ctr.SSA=2;
ctr.nsteps=101;
ctr.dt=0.01;
ctr.Tinit=0; % SET 1 FOR INITIALIZATION !
ctr.Tinv=10;
ctr.TinvMelt=0.01;
ctr.HinvMelt=10;

%KoriModel('ASE2km_initsia','ASE2km_initssa1',ctr);

%% Second SSA inversion (As, MeltInv) - long run
ctr.nsteps=100001;
ctr.dt=0.2;
ctr.TinvMelt=5;

%KoriModel('ASE2km_initssa1','ASE2km_initssa2',ctr);

%% Second SSA inversion (As, MeltInv) - with GroundedMelt
ctr.nsteps=20001;
ctr.dt=0.2;
ctr.GroundedMelt=1;
KoriModel('ASE2km_initssa2','ASE2km_initssa3',ctr);

%% Control run with optimized melt (MeltInv)
ctr.GroundedMelt=0;
ctr.inverse=0;
ctr.nsteps=1001;
ctr.meltfunc=11; % use MeltInv as melt

KoriModel('ASE2km_initssa3','ASE2km_ctrl',ctr);

%% Forcing run

ctr.calving=4;
ctr.FrontalMelt=1;
ctr.meltfunc=3;
ctr.dt=0.1;

KoriModel('ASE2km_initssa2','ASE2km_fc',ctr); % start from optimized run

ctr.damage=1;

% KoriModel('Thwaitesintc1','Thwaitesrunf1',ctr);

% ctr.GeoidCalc=1;
% ctr.BedAdj=1;
% KoriModel('Thwaitesintc','Thwaitesrune',ctr);

% ctr.damage=0;
% ctr.meltfunc=8;
% fc.butfac=zeros(ctr.nsteps,1);
% KoriModel('Thwaitesintc','Thwaitesrung',ctr,fc);


end