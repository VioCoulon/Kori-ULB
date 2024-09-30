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
save('ThuleLSF5','LSF');

%  p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 700e3);
%  XX=reshape(X, [numel(X),1]);
%  YY=reshape(X, [numel(Y),1]);
%  IceMask=inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));
%  LSF2=zeros(ctr.imax,ctr.jmax);
%  LSF2(IceMask==1)=1;
%  LSF2(IceMask==0)=-1;
%  save('ThuleLSF2','LSF2');

ctr.CF_Boundary='CircThule.mat';
ctr.LSFfile='ThuleLSF5';

H=zeros(ctr.imax,ctr.jmax)+10;
Mb=zeros(ctr.imax,ctr.jmax)+0.3;
Ts=zeros(ctr.imax,ctr.jmax)-5.0;
save('ThuleIn','B','H','Mb','Ts');

% 1. Initial spin up
ctr.shelftune=1;
ctr.SSA=1; % ctr.SSA=1
ctr.dt=2; %1 
ctr.nsteps=5000; % Jim: 10000; Daniel: 6000
ctr.timeslice=1;
ctr.snapshot=100; % Daniel:100. Jim_: nothing.
KoriModel('ThuleIn','Thule_p-t',ctr); 

% 2. Adjustment to imposition of Calving Front
%ctr.Calve_Mass=1; % Not needed in Vio's version.
%ctr.CalveGround=0; % Not needed in Vio's version.
%ctr.CalveCirc=1; % Not needed in Vio's version.
ctr.WV=0;
ctr.dt=1;  
ctr.calving=2;    % Direct, constant imposition of change in front positon.
ctr.LSFReset=30; % Jim: 50. Daniel: 30
ctr.nsteps=4000; % Jim: 6000. Daniel: 4000
%ctr.MMELT=50;

%save('Thule','LSF','-append');
%KoriModel('Thule','Thule-Exp3_Hdaniel',ctr); 

% 3. Impose zero rate of calving position change WV.
ctr.WV=0;
ctr.nsteps=100; 
ctr.timeslice=1;
%ctr.LSFReset=40;
ctr.snapshot=10; % Jim: 10; Daniel: 100
%ctr.CalveCirc=1;
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