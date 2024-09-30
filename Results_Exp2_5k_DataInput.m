clear all
close all


ExpName='CalvingMIP-Exp2-Kori.nc';

if exist(ExpName,'file')==2
    delete(ExpName);
end
%load Exp2_5_toto.mat
load Exp2_5_visceff1e10_limityes_ij_toto.mat


X=-800e3:5e3:800e3;
Y=-800e3:5e3:800e3;

X5=-802.5e3:5e3:802.5e3;
Y5=-802.5e3:5e3:802.5e3;

[Xg,Yg] = meshgrid(X5,Y5);
[Xr,Yr] = meshgrid(X,Y);
Time100=0:100:1000;


iii=0;
for t = 0:100:1000

    N='Exp2_5_visceff1e10_limityes_ij_';
    T=sprintf('%03d', t);
    NT=strcat(N,T);
    if t<1000
        load(NT)
    else
        load Exp2_5_visceff1e10_limityes_ij_tot.mat
    end
    iii=iii+1;

MASK(MASK==0)=MASK(MASK==0)+2;
MASK(LSF<0)=3;

% 
 VX(:,:,t+1)=ux;
 VY(:,:,t+1)=uy;

CR(MASK==3)=nan;
H(MASK==3)=nan;
uy(MASK==3)=nan;
ux(MASK==3)=nan;

CRsi=griddedInterpolant(Xg',Yg',CR');
Bsi=griddedInterpolant(Xg',Yg',B');
VXsi=griddedInterpolant(Xg',Yg',ux');
VYsi=griddedInterpolant(Xg',Yg',uy');
Hsi=griddedInterpolant(Xg',Yg',H');
Msi=griddedInterpolant(Xg',Yg',MASK','nearest');

VXa=VXsi(Xr',Yr');
VYa=VYsi(Xr',Yr');
Ha=Hsi(Xr',Yr');
Ma=Msi(Xr',Yr');
CRa=CRsi(Xr',Yr');

calverate(:,:,iii)=CRa;
xvelmean(:,:,iii)=VXa;
yvelmean(:,:,iii)=VYa;
lithk(:,:,iii)=Ha;
mask(:,:,iii)=Ma;
topg(:,:,iii)=Bsi(Xr',Yr');

end

Time1=0:1:1000;

iii=0;
P = readtable('Circle_Profiles.csv');

 for t = 0:1:1000

     ii=t+1;
     N='Exp2_5_';
     T=sprintf('%03d', t);
     NT=strcat(N,T);
     if t<1000
         load(NT)
     else
         load Exp2_5_visceff1e10_limityes_ij_toto.mat
     end


    N='Exp2_5_visceff1e10_limityes_ij_';
    T=sprintf('%03d', t);
    NT=strcat(N,T);
    if t<1000
        load(NT)
    else
        load Exp2_5_visceff1e10_limityes_ij_toto.mat
    end
    iii=iii+1;

MASK(MASK==0)=MASK(MASK==0)+2;
MASK(LSF<0)=3;

% 
 VX(:,:,t+1)=ux;
 VY(:,:,t+1)=uy;

CR(MASK==3)=nan;
H(MASK==3)=nan;
uy(MASK==3)=nan;
ux(MASK==3)=nan;

CRsi=griddedInterpolant(Xg',Yg',CR');
Bsi=griddedInterpolant(Xg',Yg',B');
VXsi=griddedInterpolant(Xg',Yg',ux');
VYsi=griddedInterpolant(Xg',Yg',uy');
Hsi=griddedInterpolant(Xg',Yg',H');
Msi=griddedInterpolant(Xg',Yg',MASK','nearest');

P = readtable('Circle_Profiles.csv');

lithkA(:,iii)=Hsi(P.Circle_Profile_A_X,P.Circle_Profile_A_Y);
lithkB(:,iii)=Hsi(P.Circle_Profile_B_X,P.Circle_Profile_B_Y);
lithkC(:,iii)=Hsi(P.Circle_Profile_C_X,P.Circle_Profile_C_Y);
lithkD(:,iii)=Hsi(P.Circle_Profile_D_X,P.Circle_Profile_D_Y);

xvelmeanA(:,iii)=VXsi(P.Circle_Profile_A_X,P.Circle_Profile_A_Y);
xvelmeanB(:,iii)=VXsi(P.Circle_Profile_B_X,P.Circle_Profile_B_Y);
xvelmeanC(:,iii)=VXsi(P.Circle_Profile_C_X,P.Circle_Profile_C_Y);
xvelmeanD(:,iii)=VXsi(P.Circle_Profile_D_X,P.Circle_Profile_D_Y);

yvelmeanA(:,iii)=VYsi(P.Circle_Profile_A_X,P.Circle_Profile_A_Y);
yvelmeanB(:,iii)=VYsi(P.Circle_Profile_B_X,P.Circle_Profile_B_Y);
yvelmeanC(:,iii)=VYsi(P.Circle_Profile_C_X,P.Circle_Profile_C_Y);
yvelmeanD(:,iii)=VYsi(P.Circle_Profile_D_X,P.Circle_Profile_D_Y);

maskA(:,iii)=Msi(P.Circle_Profile_A_X,P.Circle_Profile_A_Y);
maskB(:,iii)=Msi(P.Circle_Profile_B_X,P.Circle_Profile_B_Y);
maskC(:,iii)=Msi(P.Circle_Profile_C_X,P.Circle_Profile_C_Y);
maskD(:,iii)=Msi(P.Circle_Profile_D_X,P.Circle_Profile_D_Y);

topgA(:,iii)=Bsi(P.Circle_Profile_A_X,P.Circle_Profile_A_Y);
topgB(:,iii)=Bsi(P.Circle_Profile_B_X,P.Circle_Profile_B_Y);
topgC(:,iii)=Bsi(P.Circle_Profile_C_X,P.Circle_Profile_C_Y);
topgD(:,iii)=Bsi(P.Circle_Profile_D_X,P.Circle_Profile_D_Y);

sA(:,iii)=P.Circle_Profile_A_S;
sB(:,iii)=P.Circle_Profile_B_S;
sC(:,iii)=P.Circle_Profile_C_S;
sD(:,iii)=P.Circle_Profile_D_S;

lithkE(:,iii)=Hsi(P.Circle_Profile_E_X,P.Circle_Profile_E_Y);
lithkF(:,iii)=Hsi(P.Circle_Profile_F_X,P.Circle_Profile_F_Y);
lithkG(:,iii)=Hsi(P.Circle_Profile_G_X,P.Circle_Profile_G_Y);
lithkH(:,iii)=Hsi(P.Circle_Profile_H_X,P.Circle_Profile_H_Y);

xvelmeanE(:,iii)=VXsi(P.Circle_Profile_E_X,P.Circle_Profile_E_Y);
xvelmeanF(:,iii)=VXsi(P.Circle_Profile_F_X,P.Circle_Profile_F_Y);
xvelmeanG(:,iii)=VXsi(P.Circle_Profile_G_X,P.Circle_Profile_G_Y);
xvelmeanH(:,iii)=VXsi(P.Circle_Profile_H_X,P.Circle_Profile_H_Y);

yvelmeanE(:,iii)=VYsi(P.Circle_Profile_E_X,P.Circle_Profile_E_Y);
yvelmeanF(:,iii)=VYsi(P.Circle_Profile_F_X,P.Circle_Profile_F_Y);
yvelmeanG(:,iii)=VYsi(P.Circle_Profile_G_X,P.Circle_Profile_G_Y);
yvelmeanH(:,iii)=VYsi(P.Circle_Profile_H_X,P.Circle_Profile_H_Y);

maskE(:,iii)=Msi(P.Circle_Profile_E_X,P.Circle_Profile_E_Y);
maskF(:,iii)=Msi(P.Circle_Profile_F_X,P.Circle_Profile_F_Y);
maskG(:,iii)=Msi(P.Circle_Profile_G_X,P.Circle_Profile_G_Y);
maskH(:,iii)=Msi(P.Circle_Profile_H_X,P.Circle_Profile_H_Y);

topgE(:,iii)=Bsi(P.Circle_Profile_E_X,P.Circle_Profile_E_Y);
topgF(:,iii)=Bsi(P.Circle_Profile_F_X,P.Circle_Profile_F_Y);
topgG(:,iii)=Bsi(P.Circle_Profile_G_X,P.Circle_Profile_G_Y);
topgH(:,iii)=Bsi(P.Circle_Profile_H_X,P.Circle_Profile_H_Y);

sE(:,iii)=P.Circle_Profile_E_S;
sF(:,iii)=P.Circle_Profile_F_S;
sG(:,iii)=P.Circle_Profile_G_S;
sH(:,iii)=P.Circle_Profile_H_S;

 end

tendlicalvf=cfflux(1:10:end)*917;
tendligroundf=glflux(1:10:end)*917;
limnsw=IVg(1:10:end)*917;
lim=(IVg(1:10:end)+IVf(1:10:end))*917;
iareagr=Ag(1:10:end);
iareafl=Af(1:10:end);

tendlicalvf(end+1)=cfflux(end)*917;
tendligroundf(end+1)=glflux(end)*917;
limnsw(end+1)=IVg(end)*917;
lim(end+1)=(IVg(end)+IVf(end))*917;
iareagr(end+1)=Ag(end);
iareafl(end+1)=Af(end);

 X=-800e3:5e3:800e3;
 Y=-800e3:5e3:800e3;


save Exp2Kori.mat mask maskA maskB maskC maskD sA sB sC sD tendlicalvf tendligroundf limnsw lim iareagr iareafl ...
    topgA topgB topgC topgD yvelmeanA yvelmeanB yvelmeanC yvelmeanD xvelmeanA xvelmeanB xvelmeanC xvelmeanD lithk ...
    lithkA lithkB lithkC lithkD xvelmean yvelmean topg maskE maskF maskG maskH sE sF sG sH topgE topgF topgG ...
    topgH yvelmeanE yvelmeanF yvelmeanG yvelmeanH xvelmeanE xvelmeanF xvelmeanG xvelmeanH lithkE lithkF lithkG lithkH calverate












