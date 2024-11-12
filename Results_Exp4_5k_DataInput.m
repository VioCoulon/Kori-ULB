clear all
close all


ExpName='CalvingMIP-Exp4-Kori.nc';

if exist(ExpName,'file')==2
    delete(ExpName);
end
load Exp4_5_gl_toto.mat

X=-800e3:5e3:800e3;
Y=-800e3:5e3:800e3;

X5=-802.5e3:5e3:802.5e3;
Y5=-802.5e3:5e3:802.5e3;

[Xg,Yg] = meshgrid(X5,Y5);
[Xr,Yr] = meshgrid(X,Y);
Time100=0:100:1000;


iii=0;
for t = 0:100:1000

    N='Exp4_5_gl_';
    T=sprintf('%03d', t);
    NT=strcat(N,T);
    if t<1000
        load(NT)
    else
        load Exp4_5_gl_toto.mat
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
topg(:,:)=Bsi(Xr',Yr');

end

Time1=0:1:1000;

iii=0;


 for t = 0:1:1000

     ii=t+1;
     N='Exp4_5_';
     T=sprintf('%03d', t);
     NT=strcat(N,T);
     if t<1000
         load(NT)
     else
         load Exp4_5_toto.mat
     end


    N='Exp4_5_';
    T=sprintf('%03d', t);
    NT=strcat(N,T);
    if t<1000
        load(NT)
    else
        load Exp4_5_toto.mat
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

P = readtable('Caprona_Profiles.csv');

lithkCapA(:,iii)=Hsi(P.Caprona_Profile_A_X,P.Caprona_Profile_A_Y);
lithkCapB(:,iii)=Hsi(P.Caprona_Profile_B_X,P.Caprona_Profile_B_Y);
lithkCapC(:,iii)=Hsi(P.Caprona_Profile_C_X,P.Caprona_Profile_C_Y);
lithkCapD(:,iii)=Hsi(P.Caprona_Profile_D_X,P.Caprona_Profile_D_Y);

xvelmeanCapA(:,iii)=VXsi(P.Caprona_Profile_A_X,P.Caprona_Profile_A_Y);
xvelmeanCapB(:,iii)=VXsi(P.Caprona_Profile_B_X,P.Caprona_Profile_B_Y);
xvelmeanCapC(:,iii)=VXsi(P.Caprona_Profile_C_X,P.Caprona_Profile_C_Y);
xvelmeanCapD(:,iii)=VXsi(P.Caprona_Profile_D_X,P.Caprona_Profile_D_Y);

yvelmeanCapA(:,iii)=VYsi(P.Caprona_Profile_A_X,P.Caprona_Profile_A_Y);
yvelmeanCapB(:,iii)=VYsi(P.Caprona_Profile_B_X,P.Caprona_Profile_B_Y);
yvelmeanCapC(:,iii)=VYsi(P.Caprona_Profile_C_X,P.Caprona_Profile_C_Y);
yvelmeanCapD(:,iii)=VYsi(P.Caprona_Profile_D_X,P.Caprona_Profile_D_Y);

maskCapA(:,iii)=Msi(P.Caprona_Profile_A_X,P.Caprona_Profile_A_Y);
maskCapB(:,iii)=Msi(P.Caprona_Profile_B_X,P.Caprona_Profile_B_Y);
maskCapC(:,iii)=Msi(P.Caprona_Profile_C_X,P.Caprona_Profile_C_Y);
maskCapD(:,iii)=Msi(P.Caprona_Profile_D_X,P.Caprona_Profile_D_Y);

topgCapA(:,iii)=Bsi(P.Caprona_Profile_A_X,P.Caprona_Profile_A_Y);
topgCapB(:,iii)=Bsi(P.Caprona_Profile_B_X,P.Caprona_Profile_B_Y);
topgCapC(:,iii)=Bsi(P.Caprona_Profile_C_X,P.Caprona_Profile_C_Y);
topgCapD(:,iii)=Bsi(P.Caprona_Profile_D_X,P.Caprona_Profile_D_Y);

sCapA(:,iii)=P.Caprona_Profile_A_S;
sCapB(:,iii)=P.Caprona_Profile_B_S;
sCapC(:,iii)=P.Caprona_Profile_C_S;
sCapD(:,iii)=P.Caprona_Profile_D_S;

P = readtable('Halbrane_Profiles.csv');

lithkHalA(:,iii)=Hsi(P.Halbrane_Profile_A_X,P.Halbrane_Profile_A_Y);
lithkHalB(:,iii)=Hsi(P.Halbrane_Profile_B_X,P.Halbrane_Profile_B_Y);
lithkHalC(:,iii)=Hsi(P.Halbrane_Profile_C_X,P.Halbrane_Profile_C_Y);
lithkHalD(:,iii)=Hsi(P.Halbrane_Profile_D_X,P.Halbrane_Profile_D_Y);

xvelmeanHalA(:,iii)=VXsi(P.Halbrane_Profile_A_X,P.Halbrane_Profile_A_Y);
xvelmeanHalB(:,iii)=VXsi(P.Halbrane_Profile_B_X,P.Halbrane_Profile_B_Y);
xvelmeanHalC(:,iii)=VXsi(P.Halbrane_Profile_C_X,P.Halbrane_Profile_C_Y);
xvelmeanHalD(:,iii)=VXsi(P.Halbrane_Profile_D_X,P.Halbrane_Profile_D_Y);

yvelmeanHalA(:,iii)=VYsi(P.Halbrane_Profile_A_X,P.Halbrane_Profile_A_Y);
yvelmeanHalB(:,iii)=VYsi(P.Halbrane_Profile_B_X,P.Halbrane_Profile_B_Y);
yvelmeanHalC(:,iii)=VYsi(P.Halbrane_Profile_C_X,P.Halbrane_Profile_C_Y);
yvelmeanHalD(:,iii)=VYsi(P.Halbrane_Profile_D_X,P.Halbrane_Profile_D_Y);

maskHalA(:,iii)=Msi(P.Halbrane_Profile_A_X,P.Halbrane_Profile_A_Y);
maskHalB(:,iii)=Msi(P.Halbrane_Profile_B_X,P.Halbrane_Profile_B_Y);
maskHalC(:,iii)=Msi(P.Halbrane_Profile_C_X,P.Halbrane_Profile_C_Y);
maskHalD(:,iii)=Msi(P.Halbrane_Profile_D_X,P.Halbrane_Profile_D_Y);

topgHalA(:,iii)=Bsi(P.Halbrane_Profile_A_X,P.Halbrane_Profile_A_Y);
topgHalB(:,iii)=Bsi(P.Halbrane_Profile_B_X,P.Halbrane_Profile_B_Y);
topgHalC(:,iii)=Bsi(P.Halbrane_Profile_C_X,P.Halbrane_Profile_C_Y);
topgHalD(:,iii)=Bsi(P.Halbrane_Profile_D_X,P.Halbrane_Profile_D_Y);

sHalA(:,iii)=P.Halbrane_Profile_A_S;
sHalB(:,iii)=P.Halbrane_Profile_B_S;
sHalC(:,iii)=P.Halbrane_Profile_C_S;
sHalD(:,iii)=P.Halbrane_Profile_D_S;

 end

tendlicalvf=cfflux(1:20:end)*917;
tendligroundf=glflux(1:20:end)*917;
limnsw=IVg(1:20:end)*917;
lim=(IVg(1:20:end)+IVf(1:20:end))*917;
iareagr=Ag(1:20:end);
iareafl=Af(1:20:end);

tendlicalvf(end+1)=cfflux(end)*917;
tendligroundf(end+1)=glflux(end)*917;
limnsw(end+1)=IVg(end)*917;
lim(end+1)=(IVg(end)+IVf(end))*917;
iareagr(end+1)=Ag(end);
iareafl(end+1)=Af(end);

 X=-800e3:5e3:800e3;
 Y=-800e3:5e3:800e3;



save Exp4Kori.mat mask maskCapA maskCapB maskCapC maskCapD sCapA sCapB sCapC sCapD tendlicalvf tendligroundf limnsw lim iareagr iareafl ...
    topgCapA topgCapB topgCapC topgCapD yvelmeanCapA yvelmeanCapB yvelmeanCapC yvelmeanCapD xvelmeanCapA xvelmeanCapB xvelmeanCapC xvelmeanCapD lithk ...
    lithkCapA lithkCapB lithkCapC lithkCapD xvelmean yvelmean topg calverate maskHalA maskHalB maskHalC maskHalD sHalA sHalB sHalC sHalD  ...
    topgHalA topgHalB topgHalC topgHalD yvelmeanHalA yvelmeanHalB yvelmeanHalC yvelmeanHalD xvelmeanHalA xvelmeanHalB xvelmeanHalC xvelmeanHalD  ...
    lithkHalA lithkHalB lithkHalC lithkHalD











