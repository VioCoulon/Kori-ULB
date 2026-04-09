
function [betaxGL,betayGL]=GroundingLineBeta(ctr,glMASK,HAF,B,SLR,par, ...
    Ax,Ay,Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy, ...
    betax,betay)

% Analogous to GroundingLineFlux, but modifies basal friction
% instead of velocity.

betaxGL = betax;
betayGL = betay;

[jGL,iGL] = meshgrid(1:ctr.jmax,1:ctr.imax);

%% ===================== X-DIRECTION (ux > 0) =====================

q0=(Ax*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
    Asfx.^(1/ctr.m)/4^par.n);
qs=q0.^(ctr.m/(ctr.m+1));
qe0=ctr.m*(par.n+2)/(ctr.m+1);

Mgl=zeros(ctr.imax,ctr.jmax);
angnorm=zeros(ctr.imax,ctr.jmax);
nx=ones(ctr.imax,ctr.jmax);
ny=zeros(ctr.imax,ctr.jmax);
Theta=ones(ctr.imax,ctr.jmax);

glMASK1=circshift(glMASK,[0 -1]);
glMASK0=circshift(glMASK,[0 1]);
ux1=circshift(ux,[0 1]);
HAF1=circshift(HAF,[0 -1]);
B1=circshift(B,[0 -1]);
Tf1=circshift(Tf,[0 -1]);
Txx1=circshift(Txx,[0 -1]);
Tyy1=circshift(Tyy,[0 -1]);
Txy1=circshift(Txy,[0 -1]);

Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & ux>0 & ux1>0)=1;

fracgx=HAF./(HAF-HAF1);
hbgx=(1-fracgx).*B+fracgx.*B1;
hgx=(SLR-hbgx)*par.rhow/par.rho;
Mgl(hgx<0)=0;

[angnorm(Mgl==1)] = arrayfun(@(iGL,jGL) ...
    NormCalc(iGL,jGL,iGL,jGL+1,glMASK,ctr), ...
    iGL(Mgl==1),jGL(Mgl==1));

nx(Mgl==1)=cos(angnorm(Mgl==1));
ny(Mgl==1)=sin(angnorm(Mgl==1));

Theta(Mgl==1)=max(min((butfac*(Txx1(Mgl==1).*nx(Mgl==1).^2+ ...
    Tyy1(Mgl==1).*ny(Mgl==1).^2+ ...
    Txy1(Mgl==1).*nx(Mgl==1).*ny(Mgl==1)) + ...
    (1-butfac)*Tf1(Mgl==1))./Tf1(Mgl==1),1),0);

if ctr.schoof==1
    Theta=Theta.^(par.n*ctr.m/(ctr.m+1));
    ugx=qs.*hgx.^qe0.*Theta;
else
    Theta=Theta.^par.n;
    ugx=q0.*hgx.^(par.n+3).*Theta;
end

ugx(Mgl==1)=ugx(Mgl==1).*abs(nx(Mgl==1));

% ===== NEW: Compute GL beta =====
betaGL = Tf1 ./ max(ugx,1e-6);

% weighting
wx=zeros(ctr.imax,ctr.jmax);
wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-ux(Mgl==1)).*hgx(Mgl==1)./1e5));

% GL cell update
betaxGL(Mgl==1)=betax(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).* ...
    (betaGL(Mgl==1)-betax(Mgl==1));

% downstream shelf
Mgl=circshift(Mgl,[0 1]);
Mgl(glMASK1<=2)=0;
Mgl(Hmx==0)=0;

betaGL=circshift(betaGL,[0 1]);
ncx=circshift(ncorx,[0 1]);
wx=circshift(wx,[0 1]);

betaxGL(Mgl==1)=betax(Mgl==1)+ncx(Mgl==1).* ...
    (sign(wx(Mgl==1)).*(1-wx(Mgl==1)).* ...
    (betaGL(Mgl==1)-betax(Mgl==1)) + ...
    (1-sign(wx(Mgl==1))).*(0 - betax(Mgl==1)));

%% ===================== X-DIRECTION (ux < 0) =====================

% (same structure — only sign + neighbors differ)

Mgl=zeros(ctr.imax,ctr.jmax);
glMASK1=circshift(glMASK,[0 -1]);
glMASK0=circshift(glMASK,[0 -2]);
ux1=circshift(ux,[0 -1]);
HAF1=circshift(HAF,[0 -1]);
B1=circshift(B,[0 -1]);

Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & ux<0 & ux1<0)=1;

fracgx=HAF1./(HAF1-HAF);
hbgx=(1-fracgx).*B1+fracgx.*B;
hgx=(SLR-hbgx)*par.rhow/par.rho;
Mgl(hgx<0)=0;

ugx=qs.*hgx.^qe0;
betaGL = Tf ./ max(ugx,1e-6);

wx=zeros(ctr.imax,ctr.jmax);
wx(Mgl==1)=max(0,min(1,(ugx(Mgl==1)-abs(ux(Mgl==1))).*hgx(Mgl==1)./1e5));

betaxGL(Mgl==1)=betax(Mgl==1)+ncorx(Mgl==1).*wx(Mgl==1).* ...
    (betaGL(Mgl==1)-betax(Mgl==1));

%% ===================== Y-DIRECTION =====================

% Same idea applied to betay
% (structure identical to your original function)

q0=(Ay*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
    Asfy.^(1/ctr.m)/4^par.n);
qs=q0.^(ctr.m/(ctr.m+1));

Mgl=zeros(ctr.imax,ctr.jmax);
glMASK1=circshift(glMASK,[-1 0]);
glMASK0=circshift(glMASK,[1 0]);
uy1=circshift(uy,[1 0]);
HAF1=circshift(HAF,[-1 0]);
B1=circshift(B,[-1 0]);

Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & uy>0 & uy1>0)=1;

fracgy=HAF./(HAF-HAF1);
hbgy=(1-fracgy).*B+fracgy.*B1;
hgy=(SLR-hbgy)*par.rhow/par.rho;
Mgl(hgy<0)=0;

ugy=qs.*hgy.^qe0;
betaGL = Tf ./ max(ugy,1e-6);

wy=zeros(ctr.imax,ctr.jmax);
wy(Mgl==1)=max(0,min(1,(ugy(Mgl==1)-uy(Mgl==1)).*hgy(Mgl==1)./1e5));

betayGL(Mgl==1)=betay(Mgl==1)+ncory(Mgl==1).*wy(Mgl==1).* ...
    (betaGL(Mgl==1)-betay(Mgl==1));

end
