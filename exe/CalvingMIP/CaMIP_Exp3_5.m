function Thule

clear; close all;
addpath /home/daniel/models/Kori-ULB;
addpath /home/daniel/models/Kori-ULB/subroutines;

%% Initial ice sheet creation
ctr.delta = 1e3; % 5e3. Try run them at 2 km!!
ctr.imax  = 1611;  % 322 (5 km), 805 (2 km), 1611 (1 km)
ctr.jmax  = 1611;

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
ctr.Asin  = zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up
ctr.Ao    = 2.9377e-18;

%Initial LSF mask
p       = nsidedpoly(1000, 'Center', [0 0], 'Radius', 755e3);
XX      = reshape(X, [numel(X),1]);
YY      = reshape(X, [numel(Y),1]);
IceMask = inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));

LSF             = zeros(ctr.imax,ctr.jmax);
LSF(IceMask==1) = 1;
LSF(IceMask==0) = -1;
save('ThuleLSF5','LSF');


%ctr.CF_Boundary = 'CircThule.mat';
ctr.LSFfile     = 'ThuleLSF5';

H  = zeros(ctr.imax,ctr.jmax)+10;
Mb = zeros(ctr.imax,ctr.jmax)+0.3;
Ts = zeros(ctr.imax,ctr.jmax)-5.0;

%---------------------------------------
% Cut out domain along symmetry axes (quarter).
ctr.imax = (ctr.imax-1)/2+2;
ctr.jmax = (ctr.jmax-1)/2+2;
ctr.Asin = zeros(ctr.imax,ctr.jmax)+1e-7;

B   = B(ctr.imax-2:end,ctr.jmax-2:end);
H   = H(ctr.imax-2:end,ctr.jmax-2:end);
Mb  = Mb(ctr.imax-2:end,ctr.jmax-2:end);
Ts  = Ts(ctr.imax-2:end,ctr.jmax-2:end);
LSF = LSF(ctr.imax-2:end,ctr.jmax-2:end);
%---------------------------------------


save('ThuleIn','B','H','Mb','Ts');

% 1. Initial spin up.
ctr.shelftune = 1;
ctr.SSA       = 1; % ctr.SSA=1
ctr.dt        = 0.2; % 1, 2, 4. 
ctr.nsteps    = 10000; % Jim: 10000; Daniel: 6000, 5000 (dt=2)
ctr.timeslice = 1;
ctr.snapshot  = 100; % Daniel:100. Jim_: nothing.
KoriModel('ThuleIn','Thule_quarter_1km',ctr); 
%KoriModel('Thule_quarter_a','Thule_quarter_b',ctr); 

% 2. Adjustment to imposition of Calving Front.
ctr.WV       = 0;
ctr.dt       = 0.2;  % dt=1
ctr.calving  = 2;    % Direct, constant imposition of change in front positon.
ctr.LSFReset = 30; % Jim: 50. Daniel: 30
ctr.nsteps   = 1000; % Jim: 6000. Daniel: 4000
%ctr.MMELT=50;
%save('Thule_quarter_a','LSF','-append');
%KoriModel('Thule_quarter_a','Thule_quarter_a2',ctr); 

% 3. Impose zero rate of calving position change WV.
ctr.WV        = 0;
ctr.nsteps    = 100; 
ctr.timeslice = 1;
ctr.snapshot  = 10; % Jim: 10; Daniel: 100
%ctr.LSFReset=40;
%KoriModel('Thule-Exp3','Exp3_5_Hdaniel',ctr); 


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