function ThuleCirc

clear; close all;

%% Initial ice sheet creation

%CalvingMIP-Algorithim 1
ctr.delta=5e3;
ctr.imax=322; 
ctr.jmax=322;


ctr.m=3;
ctr.dt=1.0; % Jim: 0.1. Daniel: 1.0
ctr.shelf=1;
ctr.shelftune=1;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up.
ctr.Ao=2.9377e-18;
ctr.LSFfile='ThuleLSF5';
ctr.SSA=1;


ctr.calving=7;       % Jim: 5 CalveMip Periodic forcing, but it is 7 in Vio's version.
ctr.nsteps=1000;    % Daniel 500 (jsut for re-advance), 1000; Jim: 10000
ctr.timeslice=1;
ctr.LSFReset=1000;
ctr.snapshot=1000; % Jim: 1000
ctr.CR_AMP=300;      % ctr.CR_AMP is max rate of front position change
ctr.MMELT=300;
ctr.t=500;

%KoriModel('Exp1_5','Exp2_5',ctr);
%KoriModel('Exp1_5','Exp2_5',ctr);

%KoriModel('Exp1_5_pseudo-transient','Exp2_5_pseudo-transient_H-frank_calv_H-lin',ctr);
%KoriModel('Exp2_5_pseudo-transient_H-frank_retrear','Exp2_5_pseudo-transient_H-frank_advance',ctr);

% Full retreate/advance.
%KoriModel('Exp1_5_visceff1e10_limityes','Exp2_5_visceff1e10_limityes_ij',ctr);
%KoriModel('Exp1_5_visceff1e10_limityes','Exp2_5_test',ctr);

KoriModel('Exp1_5_visceff8e9','Exp2_5_visceff8e9',ctr);

% Just for re-advance starting from a fully retrated calving front.
%KoriModel('Exp2_5_visceff1e10_limityes_retreated','Exp2_5_visceff1e10_limityes_ij',ctr);

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