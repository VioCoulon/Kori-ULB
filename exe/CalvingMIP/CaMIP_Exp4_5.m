function Thule

clear; close all;

addpath /home/daniel/models/Kori-ULB;
addpath /home/daniel/models/Kori-ULB/subroutines;

parent_path = '/home/daniel/models/Kori-ULB/';
%path_in     = [parent_path, 'ice_data/calvingMIP/Exp3/dx_2km/'];
path_in     = [parent_path, 'output/CalvingMIP/Exp3/dx_5km/OceanVisc_7e9/'];
path_out    = [parent_path, 'output/CalvingMIP/Exp4/dx_5km/OceanVisc_7e9/'];


ctr.runmode  = 1;      % 1: graphics; 3: no graphics

%% Initial ice sheet creation
quarter   = true;
ctr.delta = 5e3; % 5e3. Try run them at 2 km!!
ctr.imax  = 321;  % 161 (10 km), 322 (5 km), 805 (2 km), 1611 (1 km)
ctr.jmax  = 321;

Li    = (ctr.imax-1)*ctr.delta;
Lj    = (ctr.jmax-1)*ctr.delta;
[X,Y] = meshgrid(-Lj/2:ctr.delta:Lj/2,-Li/2:ctr.delta:Li/2);

R  = 800e3 ;%800e3
Bc = 900;  %900
Bl = -2000; %-2000
Ba = 1100;  %1100
B  = BedGeom(X,Y,R,Bc,Bl,Ba);

ctr.m     = 3;
ctr.dt    = 1;
ctr.shelf = 1;
ctr.Asin  = zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmarsh set up
ctr.Ao    = 2.9377e-18;

%Initial LSF mask
p       = nsidedpoly(1000, 'Center', [0 0], 'Radius', 755e3);
XX      = reshape(X, [numel(X),1]);
YY      = reshape(X, [numel(Y),1]);
IceMask = inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));


LSF=zeros(ctr.imax,ctr.jmax);
LSF(IceMask==1)=1;
LSF(IceMask==0)=-1;
%save('ThuleLSF5','LSF');
out_1 = [path_in, 'ThuleLSF5'];
save(out_1, 'LSF');

%ctr.CF_Boundary = 'CircThule.mat';
ctr.LSFfile     = 'ThuleLSF5';

H  = zeros(ctr.imax,ctr.jmax)+10;
Mb = zeros(ctr.imax,ctr.jmax)+0.3;
Ts = zeros(ctr.imax,ctr.jmax)-5.0;

%---------------------------------------
% Cut out domain along symmetry axes (quarter).
if quarter == true

    ctr.mismip = 2;                 % Necessary for boundary conditions.
    ctr.imax   = (ctr.imax-1)/2+2;
    ctr.jmax   = (ctr.jmax-1)/2+2;
    ctr.Asin   = zeros(ctr.imax,ctr.jmax)+1e-7;

    B   = B(ctr.imax-2:end,ctr.jmax-2:end);
    H   = H(ctr.imax-2:end,ctr.jmax-2:end);
    Mb  = Mb(ctr.imax-2:end,ctr.jmax-2:end);
    Ts  = Ts(ctr.imax-2:end,ctr.jmax-2:end);
    LSF = LSF(ctr.imax-2:end,ctr.jmax-2:end);
end
%---------------------------------------


% 1. Initial spin up.
ctr.SSA=1;
%ctr.nsteps=6000;  % Jim: 10000; Daniel: 6000
%KoriModel('ThuleIn','Thule',ctr); 

% 2. Adjustment to imposition of Calving Front
ctr.MMELT=50;
ctr.Calve_Mass=1;
ctr.CalveGround=0;
ctr.WV=0;
ctr.dt=1;  
ctr.calving=2;
ctr.LSFReset=50;
ctr.nsteps=4000;
ctr.timeslice=1;
ctr.snapshot=10;
%save('Thule_p-t','LSF','-append');
%KoriModel('Thule_p-t','Thule-Circ_p-t',ctr); 

% 3. Impose zero rate of calving position change WV.
ctr.WV=0;
ctr.nsteps=100; 
ctr.timeslice=1;
ctr.LSFReset=40;
ctr.snapshot=10;
%KoriModel('Thule-Circ_p-t','Exp3_5_p-t',ctr); 


% 4. Calving MIP experiment 4 forcing.
%ctr.SSA=1;
%ctr.calving=7;
%ctr.nsteps=1000; 
%ctr.timeslice=1;
%ctr.LSFReset=50; %100, 250, 1000
%ctr.snapshot=50;    % 10000
%ctr.dt=1.0;  % Jim 0.05, Daniel: 1.0
%ctr.CR_AMP=750;

ctr.calving   = 8;
ctr.nsteps    = 1000;   % 1000 (dt=1), 10000
ctr.timeslice = 1;
ctr.LSFReset  = 50; %100, 250, 1000
ctr.snapshot  = 25;  % 1000
ctr.dt        = 1.0;  % Jim 0.05, Daniel: 1.0 (10 km), 0.1 (dx=2 km)
ctr.CR_AMP    = 750;


%KoriModel('Exp3_5_p-t','Exp4_5_gl3_visc1e7',ctr); 

out_1 = [path_in, 'Thule_quarter_Exp3_3'];
out_2 = [path_out, 'Thule_quarter_Exp4'];

KoriModel(out_1, out_2, ctr); 

%KoriModel('Exp3_5_p-t','Exp4_5_p-t_retreat',ctr); 


% START FROM RETREATED POSITION.
% Load the original LSF field (Exp3) to force limits in the readvance.
%load('Exp3_5_p-t_009','LSF')
%LSFo=LSF;
%save('Exp4_5_p-t_retreat','LSFo','-append');
%KoriModel('Exp4_5_p-t_retreat','Exp4_5_p-t_advanced',ctr); 


end

function [B]=BedGeom(x,y,R,Bc,Bl,Ba)
% param ters

rc=0;
%polarcoordinates
r=sqrt(x.*x+y.*y);
theta=atan2(y,x);
% B calculation
l=R-cos(2*theta).*R/2;
a=Bc-(Bc-Bl)*(r-rc).^2./(R-rc).^2;
B=Ba*cos(3*pi*r./l)+a ;
end