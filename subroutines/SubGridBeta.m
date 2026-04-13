
function [betaxGL,betayGL]=SubGridBeta(ctr,glMASK,HAF,B,SLR,par, ...
    Ax,Ay,Tf,Txx,Tyy,Txy,butfac,ncorx,ncory,ux,uy,Hmx,Hmy,Asfx,Asfy, ...
    betax,betay,cnt)

% Analogous to GroundingLineFlux, but modifies basal friction
% instead of velocity.

betaxGL = betax;
betayGL = betay;

uxsch=ux;
uysch=uy;

[jGL,iGL] = meshgrid(1:ctr.jmax,1:ctr.imax); %VL: grid for GL orientation

% GL flux in x-direction
q0  = (Ax*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
        Asfx.^(1/ctr.m)/4^par.n);
qs  = q0.^(ctr.m/(ctr.m+1));
qe0 = ctr.m*(par.n+2)/(ctr.m+1);

% ------------------------------------------------------------
% ------------------------------------------------------------
% conditions for GL in x-direction and ux>0
qgx     = zeros(ctr.imax,ctr.jmax);
Mgl     = zeros(ctr.imax,ctr.jmax);
angnorm = zeros(ctr.imax,ctr.jmax);   %VL
nx      = ones(ctr.imax,ctr.jmax); %VL
ny      = zeros(ctr.imax,ctr.jmax); %VL
Theta   = ones(ctr.imax,ctr.jmax); %VL

glMASK1 = circshift(glMASK,[0 -1]); % glMASK(i,j+1)
glMASK0 = circshift(glMASK,[0 1]); % glMASK(i,j-1), upstream
ux1     = circshift(ux,[0 1]); % ux(i,j-1), upstream velocity
betax1  = circshift(betax,[0 1]); % betax(i,j-1), upstream beta
HAF1    = circshift(HAF,[0 -1]); % HAF(i,j+1)
B1      = circshift(B,[0 -1]); % B(i,j+1)
Tf1     = circshift(Tf,[0 -1]); %VL:  Tf(i,j+1)
Txx1    = circshift(Txx,[0 -1]); %VL:  Txx(i,j+1)
Tyy1    = circshift(Tyy,[0 -1]); %VL:  Tyy(i,j+1)
Txy1    = circshift(Txy,[0 -1]); %VL:  Txy(i,j+1)

Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & ux>0 & ux1>0) = 1; % u-grid(i,j) 

fracgx = HAF./(HAF-HAF1);
hbgx   = (1.-fracgx).*B+fracgx.*B1;
hgx    = (SLR-hbgx)*par.rhow/par.rho;

Mgl(hgx<0) = 0;

% Handy definiton for mask.
a = Mgl==1;

% GL in x-direction, u>0 --> jg=jGL, jh=jg+1, ih=ig=iGL
[angnorm(a)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL,jGL+1, ...
                                glMASK,ctr),iGL(a),jGL(a));

nx(a) = cos(angnorm(a));
ny(a) = sin(angnorm(a));

Theta(a) = max(min((butfac*(Txx1(a).*nx(a).^2+Tyy1(a).* ...
                            ny(a).^2+Txy1(a).*nx(a).*ny(a))+(1.-butfac)* ...
                                Tf1(a))./Tf1(a),1),0);

%if ctr.schoof==1 % Schoof%
%    Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
%else % Pattyn
%    Theta=Theta.^(par.n);
%end
%if ctr.schoof==1 % Schoof
%    ugx=qs.*hgx.^qe0.*Theta;
%else % Pattyn
%    ugx=q0.*hgx.^(par.n+3).*Theta;
%end

% Daniel.
Theta = Theta.^(par.n*ctr.m/(ctr.m+1.));
ugx   = qs.*hgx.^qe0.*Theta;

ugx(a) = ugx(a).*abs(nx(a));
qgx(a) = ugx(a).*hgx(a);

% weighting factor (Pollard & DeConto, 2020)
wx    = zeros(ctr.imax,ctr.jmax); 
wx(a) = max(0,min(1,(ugx(a)-ux(a)).*hgx(a)./10^5));

% GL grid cell
uxsch(a) = ux(a)+ncorx(a).*wx(a).*(qgx(a)./max(Hmx(a),0.01)-ux(a));

% ------------------------------------------------------------
% UPDATE BETA WITH SCHOOF VELOCITIES.
if ctr.uSSAexist==1 || cnt>1
    ussa=vec2h(uxsch,uysch);    %VL: ussa on h-grid
else
    ussa=zeros(ctr.imax,ctr.jmax)+0.1;
end
ussa=max(ussa,1e-3);
if par.ShelfPinning==1 && ctr.stdBexist==1 && ctr.inverse==0 
    fg=max(0,1-(HB-B)./stdB); % subgrid pinning points in ice shelf
else
    fg=1;
end
if ctr.u0>1e10
    beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
else
    beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).* ...
        Asf/ctr.u0).^(-1/ctr.m);
end
if ctr.SSA==3
    % Integration factor.
    [F2,~,~]=Fint(ctr,etaD,H,zeta);
    % Effective beta from F2 correction in DIVA model.
    beta2=beta2./(1.0+beta2.*F2);
    beta2(beta2>1e8)=1./F2(beta2>1e8); % frozen conditions
end
%beta2=min(beta2,1e8);
%beta2(MASK==0)=0;
betagx=0.5*(beta2+circshift(beta2,[0 -1]));
%betagy=0.5*(beta2+circshift(beta2,[-1 0]));
% ------------------------------------------------------------


% Daniel correction on beta (at the grounding line).
betagx(a)   = betagx(a).*abs(nx(a)).*hgx(a) ./ max(Hmx(a),0.01);
betaxsch(a) = betax(a) + ncorx(a).*wx(a).*(betagx(a)-betax(a));

                    
% Ice shelf grid cell downstream of GL
Mgl=circshift(Mgl,[0 1]); % condition for i,j+1
Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
Mgl(Hmx==0)=0; % only apply if shelf exists

% Handy definiton for mask.
a = Mgl==1;

qgx = circshift(qgx,[0 1]);
ncx = circshift(ncorx,[0 1]);
wx  = circshift(wx,[0 1]);

betagx(a)   = betagx(a).*abs(nx(a)).*hgx(a) ./ max(Hmx(a),0.01);

uxsch(a)=ux(a)+ncx(a).*( ...
                sign(wx(a)).*(1-wx(a)).*(qgx(a)./ ...
                    max(Hmx(a),0.01)-ux(a)) ... % wx>0 (retreat)
                        +(1-sign(wx(a))).*(min(qgx(a)./ ...
                            max(Hmx(a),0.01),ux1(a))-ux(a)));   % wx=0 (advance)

%betax(a)=betax(a)+ncx(a).*( ...
%                sign(wx(a)).*(1-wx(a)).*(betagx(a)-betax(a)) ... % wx>0 (retreat)
%                        +(1-sign(wx(a))).*(min(betagx(a),betax1(a))-betax(a)));   % wx=0 (advance)


% ------------------------------------------------------------
% ------------------------------------------------------------



% ------------------------------------------------------------
% ------------------------------------------------------------
% conditions for GL in x-direction and ux<0
qgx     = zeros(ctr.imax,ctr.jmax);
Mgl     = zeros(ctr.imax,ctr.jmax);
angnorm = zeros(ctr.imax,ctr.jmax);   %VL
nx      = ones(ctr.imax,ctr.jmax); %VL
ny      = zeros(ctr.imax,ctr.jmax); %VL
Theta   = ones(ctr.imax,ctr.jmax); %VL

glMASK1 = circshift(glMASK,[0 -1]); % glMASK(i,j+1)
glMASK0 = circshift(glMASK,[0 -2]); % glMASK(i,j+2), upstream
ux1     = circshift(ux,[0 -1]);   % ux(i,j+1), upstream velocity
betax1  = circshift(betax,[0 -1]);   % betax(i,j+1), upstream betas
HAF1    = circshift(HAF,[0 -1]); % HAF(i,j+1)
B1      = circshift(B,[0 -1]); % B(i,j+1)

Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & ux<0 & ux1<0)=1; % u-grid(i,j)

fracgx = HAF1./(HAF1-HAF);
hbgx   = (1.-fracgx).*B1+fracgx.*B;
hgx    = (SLR-hbgx)*par.rhow/par.rho;
Mgl(hgx<0)=0;

% Handy definiton for mask.
a = Mgl==1;

%VL: GL in x-direction, u<0 --> jh=jGL, jg=jh+1, ih=ig=iGL
[angnorm(a)] = arrayfun(@(iGL,jGL) NormCalc(iGL,jGL+1,iGL,jGL, ...
                glMASK,ctr),iGL(a),jGL(a));  %VL

nx(a) = cos(angnorm(a));
ny(a) = sin(angnorm(a));

Theta(a) = max(min((butfac*(Txx(a).*nx(a).^2+Tyy(a).* ...
                ny(a).^2+Txy(a).*nx(a).*ny(a))+(1.-butfac)* ...
                    Tf(a))./Tf(a),1),0);

%if ctr.schoof==1 % Schoof
%    Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
%else % Pattyn
%    Theta=Theta.^(par.n);
%end
%if ctr.schoof==1 % Schoof
%    ugx=qs.*hgx.^qe0.*Theta;
%else % Pattyn
%    ugx=q0.*hgx.^(par.n+3).*Theta;
%end

% Daniel.
Theta = Theta.^(par.n*ctr.m/(ctr.m+1.));
ugx   = qs.*hgx.^qe0.*Theta;

ugx(a) = ugx(a).*abs(nx(a));
qgx(a) = -ugx(a).*hgx(a);

% weighting factor (Pollard & DeConto, 2020)
wx    = zeros(ctr.imax,ctr.jmax); 
wx(a) = max(0,min(1,(ugx(a)-abs(ux(a))).*hgx(a)./10^5));

% GL grid cell
uxsch(a)=ux(a)+ncorx(a).*wx(a).*(qgx(a)./max(Hmx(a),0.01)-ux(a));




% ice shelf grid cell downstream of GL
Mgl=circshift(Mgl,[0 -1]); % condition for i,j-1
Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
Mgl(Hmx==0)=0; % only apply if shelf exists

% Handy definiton for mask.
a = Mgl==1;

qgx=circshift(qgx,[0 -1]);
ncx=circshift(ncorx,[0 -1]);
wx=circshift(wx,[0 -1]);

uxsch(a)=ux(a)+ncx(a).*( ...
    sign(wx(a)).*(1-wx(a)).*(qgx(a)./ ...
    max(Hmx(a),0.01)-ux(a)) ... % wx>0 (retreat)
    +(1-sign(wx(a))).*(max(qgx(a)./ ...
    max(Hmx(a),0.01),ux1(a))-ux(a)));   % wx=0 (advance)

% ------------------------------------------------------------
% ------------------------------------------------------------



% GL flux in y-direction
q0=(Ay*(par.rho*par.g)^(par.n+1)*(1-par.rho/par.rhow)^par.n.* ...
    Asfy.^(1/ctr.m)/4^par.n);
qs=q0.^(ctr.m/(ctr.m+1));

% conditions for GL in y-direction and uy>0
qgx     = zeros(ctr.imax,ctr.jmax);
Mgl     = zeros(ctr.imax,ctr.jmax);
angnorm = zeros(ctr.imax,ctr.jmax);   %VL
nx      = ones(ctr.imax,ctr.jmax); %VL
ny      = zeros(ctr.imax,ctr.jmax); %VL
Theta   = ones(ctr.imax,ctr.jmax); %VL

glMASK1 = circshift(glMASK,[-1 0]); % glMASK(i+1,j)
glMASK0 = circshift(glMASK,[1 0]); % glMASK(i-1,j), upstream
uy1     = circshift(uy,[1 0]);    % uy(i-1,j), upstream velocity
HAF1    = circshift(HAF,[-1 0]); % HAF(i+1,j)
B1      = circshift(B,[-1 0]); % B(i+1,j)
Tf1     = circshift(Tf,[-1 0]); %VL:  Tf(i+1,j)
Txx1    = circshift(Txx,[-1 0]); %VL:  Txx(i+1,j)
Tyy1    = circshift(Tyy,[-1 0]); %VL:  Tyy(i+1,j)
Txy1    = circshift(Txy,[-1 0]); %VL:  Txy(i+1,j)

Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & uy>0 & uy1>0)=1; % v-grid(i,j)
fracgy=HAF./(HAF-HAF1);
hbgy=(1.-fracgy).*B+fracgy.*B1;
hgy=(SLR-hbgy)*par.rhow/par.rho;
Mgl(hgy<0)=0;

% Handy definiton for mask.
a = Mgl==1;

%VL: GL in y-direction, v>0 --> ig=iGL, ih=ig+1, jh=jg=jGL
[angnorm(a)]=arrayfun(@(iGL,jGL) NormCalc(iGL,jGL,iGL+1,jGL, ...
    glMASK,ctr),iGL(a),jGL(a));  %VL
nx(a)=cos(angnorm(a));
ny(a)=sin(angnorm(a));
Theta(a)=max(min((butfac*(Txx1(a).*nx(a).^2+ ...
    Tyy1(a).*ny(a).^2+Txy1(a).*nx(a).*ny(a))+ ...
    (1.-butfac)*Tf1(a))./Tf1(a),1),0);
if ctr.schoof==1 % Schoof
    Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
else % Pattyn
    Theta=Theta.^(par.n);
end
if ctr.schoof==1 % Schoof
    ugy=qs.*hgy.^qe0.*Theta;
else % Pattyn
    ugy=q0.*hgy.^(par.n+3).*Theta;
end
ugy(a)=ugy(a).*abs(ny(a));
qgy(a)=ugy(a).*hgy(a);
wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
wy(a)=max(0,min(1,(ugy(a)-uy(a)).*hgy(a)./10^5));

% GL grid cell
uysch(a)=uy(a)+ncory(a).*wy(a).*(qgy(a)./ ...
    max(Hmy(a),0.01)-uy(a));

    % ice shelf grid cell downstream of GL
Mgl=circshift(Mgl,[1 0]); % condition for i+1,j
Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
Mgl(Hmy==0)=0; % only apply if shelf exists

% Handy definiton for mask.
a = Mgl==1;


qgy=circshift(qgy,[1 0]);
ncy=circshift(ncory,[1 0]);
wy=circshift(wy,[1 0]);
uysch(a)=uy(a)+ncy(a).*( ...
    sign(wy(a)).*(1-wy(a)).*(qgy(a)./ ...
    max(Hmy(a),0.01)-uy(a)) ...   % wy>0 (retreat)
    +(1-sign(wy(a))).*(min(qgy(a)./ ...
    max(Hmy(a),0.01),uy1(a))-uy(a)));   % wy=0 (advance)

% conditions for GL in y-direction and uy<0
qgx     = zeros(ctr.imax,ctr.jmax);
Mgl     = zeros(ctr.imax,ctr.jmax);
angnorm = zeros(ctr.imax,ctr.jmax);   %VL
nx      = ones(ctr.imax,ctr.jmax); %VL
ny      = zeros(ctr.imax,ctr.jmax); %VL
Theta   = ones(ctr.imax,ctr.jmax); %VL

glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
glMASK0=circshift(glMASK,[-2 0]); % glMASK(i+2,j), upstream
uy1=circshift(uy,[-1 0]);   % uy(i+1,j), upstream velocity
HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
B1=circshift(B,[-1 0]); % B(i+1,j)

Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & uy<0 & uy1<0)=1; % v-grid(i,j)
fracgy=HAF1./(HAF1-HAF);
hbgy=(1.-fracgy).*B1+fracgy.*B;
hgy=(SLR-hbgy)*par.rhow/par.rho;
Mgl(hgy<0)=0;

% Handy definiton for mask.
a = Mgl==1;


%VL: GL in y-direction, v<0 --> ih=iGL, ig=ih+1, jh=jg=jGL
[angnorm(a)]=arrayfun(@(iGL,jGL) NormCalc(iGL+1,jGL,iGL,jGL, ...
    glMASK,ctr),iGL(a),jGL(a));  %VL
nx(a)=cos(angnorm(a));
ny(a)=sin(angnorm(a));
Theta(a)=max(min((butfac*(Txx(a).*nx(a).^2+Tyy(a).* ...
    ny(a).^2+Txy(a).*nx(a).*ny(a))+(1.-butfac)* ...
    Tf(a))./Tf(a),1),0);
if ctr.schoof==1 % Schoof
    Theta=Theta.^(par.n*ctr.m/(ctr.m+1.));
else % Pattyn
    Theta=Theta.^(par.n);
end
if ctr.schoof==1 % Schoof
    ugy=qs.*hgy.^qe0.*Theta;
else % Pattyn
    ugy=q0.*hgy.^(par.n+3).*Theta;
end

ugy(a)=ugy(a).*abs(ny(a));
qgy(a)=-ugy(a).*hgy(a);
wy=zeros(ctr.imax,ctr.jmax); % weighting factor (Pollard & DeConto, 2020)
wy(a)=max(0,min(1,(ugy(a)-abs(uy(a))).*hgy(a)./10^5));
% GL grid cell
uysch(a)=uy(a)+ncory(a).*wy(a).*(qgy(a)./ ...
    max(Hmy(a),0.01)-uy(a));

    % ice shelf grid cell downstream of GL
Mgl=circshift(Mgl,[-1 0]); % condition for i-1,j
Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
Mgl(Hmy==0)=0; % only apply if shelf exists

% Handy definiton for mask.
a = Mgl==1;

qgy=circshift(qgy,[-1 0]);
ncy=circshift(ncory,[-1 0]);
wy=circshift(wy,[-1 0]);

uysch(a)=uy(a)+ncy(a).*( ...
    sign(wy(a)).*(1-wy(a)).*(qgy(a)./ ...
    max(Hmy(a),0.01)-uy(a)) ...   % wy>0 (retreat)
    +(1-sign(wy(a))).*(max(qgy(a)./ ...
    max(Hmy(a),0.01),uy1(a))-uy(a)));   % wy=0 (advance)

%cnt


end
