function [db]=BasalDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,HAF)
    
    % Kori-ULB
    % Basal Damage functions

    % ctr.bsldamage=0: No damage
    % ctr.bsldamage=1: Nye damage function      (following Sun et al., 2017: 10.5194/tc-11-2543-2017)
    % ctr.bsldamage=2: Weertman damage function (following Lai et al., 2020: 10.1038/s41586-020-2627-8)
    %                  Is the factor pi/2 valid for basal crevasses?
    % ctr.bsldamage=3: Kachuck damage function  (following Kachuck et al., 2022: 10.1017/jog.2022.12 )

    eps=1e-8; % avoid zero values
    [lambda1,lambda2]=PrincipalStrain(dudx,dvdy,dudy,dvdx); % 1st/2nd principal strain
    % convert strain to stress -- note that eta is the vertically integrated viscosity
    % hence, the ice thickness is considered in there
    tau1=2*lambda1.*eta./(H+eps);

    if ctr.bsldamage~=0
        if ctr.bsldamage==1
            db=(par.rho/(par.rhow-par.rho))*((tau1./(par.rho*par.g))-max(HAF,0));
        elseif ctr.bsldamage==2
            db=(par.rho/(par.rhow-par.rho))*((pi*0.5*tau1./(par.rho*par.g))-max(HAF,0));
        elseif ctr.bsldamage==3
            alpha=lambda2./lambda1;
            db=db.*(2+alpha);
            %    db=(par.rho/(par.rhow-par.rho))*(tau1.*(2+alpha)./(par.rho*par.g*(H+eps)));
        end
    else
        db=zeros(ctr.imax,ctr.jmax); % Initialize to zeros
    end

    % CWR: crevasse width ratio
    % CW:  crevasse width [m]
    CWR = ctr.CrWidth/ctr.delta;
    db=max(0.0,db.*CWR);
    db(tau1<ctr.tauice)=0; % no damage if yield strength of ice is not reached

    % db cannot be bigger than ice thickness
    db=max(0,min(db,H));

end
