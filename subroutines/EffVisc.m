function [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,MASK, ...
    glMASK,shelftune,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined based on a fixed ice thickness (Hshelf) that
% determines eta1

    dudx=(uxssa-circshift(uxssa,[0 1]))/ctr.delta;
    % Daniel choice: SYMMETRIC!!
    %dudx=0.5*(circshift(uxssa,[0 -1])-circshift(uxssa,[0 1]))/ctr.delta;
    dudx(:,1)=dudx(:,2);
    dudx(:,ctr.jmax)=dudx(:,ctr.jmax-1);
    dvdy=(uyssa-circshift(uyssa,[1 0]))/ctr.delta;
    % Daniel
    %dvdy=0.5*(circshift(uyssa,[-1 0])-circshift(uyssa,[1 0]))/ctr.delta;
    dvdy(1,:)=dvdy(2,:);
    dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
    dudy=(circshift(uxssa,[-1 1])+circshift(uxssa,[-1 0])- ...
        circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(4*ctr.delta);
    % Daniel.
    %dudy=(uxssa+circshift(uxssa,[0 1])- ...
    %    circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(2*ctr.delta);
    dudy(1,:)=dudy(2,:);
    dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
    dudy(:,1)=dudy(:,2);
    dudy(:,ctr.jmax)=dudy(:,ctr.jmax-1);
    dvdx=(circshift(uyssa,[0 -1])+circshift(uyssa,[1 -1])- ...
        circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(4*ctr.delta);
    % Daniel.
    %dvdx=(uyssa+circshift(uyssa,[1 0])- ...
    %    circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(2*ctr.delta);
    dvdx(1,:)=dvdx(2,:);
    dvdx(ctr.imax,:)=dvdx(ctr.imax-1,:);
    dvdx(:,1)=dvdx(:,2);
    dvdx(:,ctr.jmax)=dvdx(:,ctr.jmax-1);
    
    if ctr.mismip>=1
        dvdx(:,1)=-dvdx(:,2);
        dvdx(:,ctr.jmax)=dvdx(:,ctr.jmax-1);
        dvdy(1,:)=dvdy(3,:);
        dudy(1,:)=dudy(3,:);
        if ctr.mismip==1
            dvdy(ctr.imax,:)=dvdy(ctr.imax-3,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-3,:);
        else
            dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
            dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
        end
    end
    EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    EffStr=max(EffStr,1e-12);
    eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune

    % Jablasco regularizitaion implementation.
    % Bassis et al., (2021) regularization
    if ctr.bassis_reg==1
        % Glen flow law
        eta_glen=eta;
        % Diffussion creep
        d_grain=5e-3; % tunable parameter? Bassis found low effect.
        eta_diff=H.*(A.^(-1/par.n))./(2*d_grain^2);
        % Plastic regime
        eta_plas = H.*(ctr.tauice)./(2*EffStr.^0.5);
        % minimum viscosity necessary for numerical convergence
        eta_min=1e10; % tunable parameter?
        % regularized viscosity
        eta=(eta_min+((eta_diff.^-1) + (eta_glen.^-1) + (eta_plas.^-1)).^-1);
    end
    
    if ctr.shelf==1 || ctr.schoof>0
        MASKb=ones(ctr.imax,ctr.jmax); % use constant eta on edges of ice shelf
        if ctr.mismip==0
            MASKb(3:ctr.imax-2,3:ctr.jmax-2)=0;
        elseif ctr.mismip==1
            MASKb(:,1:ctr.jmax-2)=0;
        elseif ctr.mismip==2
            MASKb(1:ctr.imax-2,1:ctr.jmax-2)=0;
        end
        MASKb(MASK==1)=0;
        % remove outliers (10%) - especially important for basins where
        % grounded parts may exist on the domain boundary
        eta(MASKb==1)=trimmean(eta(MASKb==1),10);
        
        % Instead of calculating effective viscosity on the sea ice,
        % keep constant viscosity. Need to further check how to deal
        % with this. May be quoted
        eta(glMASK==6)=8.0e9; % Default: 1.0e7. Daniel: 1.0e10
    end
    
end


