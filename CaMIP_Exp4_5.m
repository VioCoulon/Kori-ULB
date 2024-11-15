function Thule

clear; close all;

%% Initial ice sheet creation

ctr.delta=5e3;
ctr.imax=322; 
ctr.jmax=322;

Li=(ctr.imax-1)*ctr.delta;
Lj=(ctr.jmax-1)*ctr.delta;
[X,Y]=meshgrid(-Lj/2:ctr.delta:Lj/2,-Li/2:ctr.delta:Li/2);

R=800e3 ;%800e3
Bc=900;  %900
Bl=-2000; %-2000
Ba=1100;  %1100
B=BedGeom(X,Y,R,Bc,Bl,Ba);

ctr.m=3;
ctr.dt=1;
ctr.shelf=1;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up
ctr.Ao=2.9377e-18;
%Initial LSF mask
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 755e3);
XX=reshape(X, [numel(X),1]);
YY=reshape(X, [numel(Y),1]);
IceMask=inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));
LSF=zeros(ctr.imax,ctr.jmax);
LSF(IceMask==1)=1;
LSF(IceMask==0)=-1;
save('ThuleLSF','LSF');

%  p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 700e3);
%  XX=reshape(X, [numel(X),1]);
%  YY=reshape(X, [numel(Y),1]);
%  IceMask=inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));
%  LSF2=zeros(ctr.imax,ctr.jmax);
%  LSF2(IceMask==1)=1;
%  LSF2(IceMask==0)=-1;
%  save('ThuleLSF2','LSF2');

%ctr.CF_Boundary='CircThule.mat';
ctr.LSFfile='ThuleLSF5';

H=zeros(ctr.imax,ctr.jmax)+10;
Mb=zeros(ctr.imax,ctr.jmax)+0.3;
Ts=zeros(ctr.imax,ctr.jmax)-5.0;
save('ThuleIn','B','H','Mb','Ts');

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
ctr.SSA=1;
ctr.calving=7;
ctr.nsteps=1000; 
ctr.timeslice=1;
ctr.LSFReset=50; %100, 250, 1000
ctr.snapshot=1000;
ctr.dt=1.0;  % Jim 0.05, Daniel: 1.0
ctr.CR_AMP=750;

KoriModel('Exp3_5_p-t','Exp4_5_gl3_visc1e7',ctr); 

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