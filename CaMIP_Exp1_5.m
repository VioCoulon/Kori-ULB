function CaMIP_Exp1_5
    
    clear; close all;
    
    %% Initial ice sheet creation.
    ctr.delta=40e3; % 5e3, 20e3
    ctr.imax=41; % 322, 82
    ctr.jmax=41; % 322, 32
    
    
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
    ctr.shelftune=1;
    ctr.SSA=1;
    ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up
    ctr.Ao=2.9377e-18;
    
    % Initial LSF mask. Prepare other fields to run Kori.
    p = nsidedpoly(10000, 'Center', [0 0], 'Radius', 750e3);
    XX=reshape(X, [numel(X),1]);
    YY=reshape(X, [numel(Y),1]);
    IceMask=inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));
    LSF=zeros(ctr.imax,ctr.jmax);
    LSF(IceMask==1)=1;
    LSF(IceMask==0)=-1;
    save('ThuleLSF5','LSF');
    H=zeros(ctr.imax,ctr.jmax)+10;
    Mb=zeros(ctr.imax,ctr.jmax)+0.3;
    Ts=zeros(ctr.imax,ctr.jmax)-5.0;
    save('ThuleIn5','B','H','Mb','Ts','LSF');
    
    % 1 - Initial spin up
    ctr.dt=1;   % Jim: 1. Daniel: 5
    ctr.nsteps=8000; % Jim: 15000, Daniel: 8000 is enough
    %KoriModel('ThuleIn5','Thule5_visceff1e10_limitno_daniel',ctr);
    ctr.timeslice=1;
    ctr.snapshot=50; % 800 output every 10 years.
    ctr.ItSolv=0; % makes a difference
    ctr.upstream=0; % makes a big difference: 0 more precise
    %KoriModel('ThuleIn5','Thule5_u-Frank_H-daniel',ctr);
    KoriModel('ThuleIn5','Thule5_upstream-0',ctr);
    
    % 2 - Adjustment to imposition of Calving Front
    ctr.calving=2;   % Direct, constant imposition of change in front positon.
    ctr.WV=0;        % ctr.WV=0 will fix calving front position to be unmoving.
    ctr.dt=1;   % 1
    ctr.LSFReset=30;
    ctr.nsteps=4000; %Jim: 10000, Daniel: 4000 is enough
    ctr.snapshot=50;
    %save('Thule5_u-Frank_H-daniel','LSF','-append'); % save('Thule5','LSF','-append'); % Make sure to save the LSF that comes from the initial spinup. 
    %KoriModel('Thule5_u-Frank_H-daniel','Thule-Circ5_pseudo-transient',ctr);
    
    % Old Jim's code.
    % %CalvingMIP-Algorithim 1
    %
    ctr.WV=0;
    ctr.nsteps=100;
    ctr.timeslice=1;
    ctr.snapshot=10;
    %KoriModel('Thule-Circ5_pseudo-transient','Exp1_5_pseudo-transient',ctr); %KoriModel('Thule-Circ5','Exp1_5',ctr);
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
B=a ;
end