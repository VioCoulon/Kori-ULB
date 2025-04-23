function [ds]=SurfaceDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,MASK)

% Kori-ULB
% Surface Damage functions

% ctr.srfdamage=0: No damage
% ctr.srfdamage=1: Nye damage function      (following Sun et al., 2017, Nick et al., 2010)
% ctr.srfdamage=2: Weertman damage function (following Lai et al., 2020)
% ctr.srfdamage=3: Kachuck based damage (removing floating contr.)
% ctr.srfdamage=4: Lai damage function      (following Lai et al., 2020)

eps=1e-8; % avoid zero values
[lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
% convert strain to stress -- note that eta is the vertically integrated viscosity
% hence, the ice thickness is considered in there
tau1=2*lambda1.*eta./(H+eps);

if ctr.srfdamage~=0

    dw=zeros(ctr.imax,ctr.jmax); % Water depth in the surface crevasse (Sun2017, Nick2010) -- TO DO!

    if ctr.srfdamage==1 % 
        ds=tau1./(par.rho*par.g)+par.rhow.*dw/par.rho;
    elseif ctr.srfdamage==2
        ds=pi*0.5*tau1./(par.rho*par.g)+par.rhow*dw/par.rho;
    elseif ctr.srfdamage==3
        alpha=lambda2./lambda1;
        ds=ds.*(2+alpha);
        %    ds=tau1.*(2+alpha)./(par.rho*par.g*(H+eps))+par.rhow*dw/par.rho;
    elseif ctr.srfdamage==4
        F=1.122;
        f=1.068;
        ds=(tau1.*pi*F)./(6*par.rho*par.g*f);
    end

else
    ds=zeros(ctr.imax,ctr.jmax); % Initialize to zeros
end

% CWR: crevasse width ratio
% CW:  crevasse width [m]
CWR = ctr.CrWidth/ctr.delta;
ds=max(0.0,ds.*CWR);
ds(tau1<ctr.tauice)=0; % no damage if yield strength of ice is not reached

% ds cannot be bigger than ice thickness
ds=max(0,min(ds,H));

end
