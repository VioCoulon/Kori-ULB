function [eta,dudx,dvdy,dudy,dvdx,d_grain,EffStr,eta_diva,dudz,dvdz]=EffViscDIVA(A,uxssa,uyssa,H,par,MASK, ...
    glMASK,shelftune,zeta,tmp,betax,betay,eta_diva,cnt,uxb_diva,uyb_diva,ctr)

% Kori-ULB
% Effective viscosity of the SSA solution. On the borders of the domain, a
% fixed value is determined based on a fixed ice thickness (Hshelf) that
% determines eta1

    dudx=(uxssa-circshift(uxssa,[0 1]))/ctr.delta;
    dudx(:,1)=dudx(:,2);
    dudx(:,ctr.jmax)=dudx(:,ctr.jmax-1);
    dvdy=(uyssa-circshift(uyssa,[1 0]))/ctr.delta;
    dvdy(1,:)=dvdy(2,:);
    dvdy(ctr.imax,:)=dvdy(ctr.imax-1,:);
    dudy=(circshift(uxssa,[-1 1])+circshift(uxssa,[-1 0])- ...
        circshift(uxssa,[1 1])-circshift(uxssa,[1 0]))/(4*ctr.delta);
    
    dudy(1,:)=dudy(2,:);
    dudy(ctr.imax,:)=dudy(ctr.imax-1,:);
    dudy(:,1)=dudy(:,2);
    dudy(:,ctr.jmax)=dudy(:,ctr.jmax-1);
    dvdx=(circshift(uyssa,[0 -1])+circshift(uyssa,[1 -1])- ...
        circshift(uyssa,[0 1])-circshift(uyssa,[1 1]))/(4*ctr.delta);

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
    % Original code for SSA.
    %EffStr=dudx.^2+dvdy.^2+dudx.*dvdy+0.25*(dudy+dvdx).^2;
    %EffStr=max(EffStr,1e-12);
    %eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    %eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Daniel: DIVA velocity implementation based on Lipscomb et al. (2019).
    % Differences in vertical dimension.
    % 3D strain rates from DIVA approximation.

    % Basal stress to compute full 3D viscosity eta(x,y,z).
    % betax and betax are defined on the velocity grid.
    %taux = zeros([ctr.imax, ctr.jmax]);
    %tauy = zeros([ctr.imax, ctr.jmax]);

    %size(uxb_diva)
    %taux = betax .* uxssa;
    %tauy = betay .* uyssa;
    taux = betax .* uxb_diva;
    tauy = betay .* uyb_diva;
    
    %size(taux)
    % Stagger tau values to the ice thickness grid (eta) to divide by eta_diva.
    taux = 0.5 * ( taux + circshift(taux, [0 1]) );
    tauy = 0.5 * ( tauy + circshift(tauy, [1 0]) );

    %any(taux(:))

    % Eq. 36 Lipscomb et al. (2019).
    %dudz = zeros([ctr.imax, ctr.jmax, ctr.kmax]);
    %dvdz = zeros([ctr.imax, ctr.jmax, ctr.kmax]);
    %zeta_m = zeros([ctr.imax, ctr.jmax, ctr.kmax]);

    % By definition: zeta(1)=0, zeta(ctr.kmax)=1;
    % Vectorial form.
    H_diva    = repmat(H, [1, 1, ctr.kmax]);
    taux_diva = repmat(taux, [1, 1, ctr.kmax]);
    tauy_diva = repmat(tauy, [1, 1, ctr.kmax]);

    k_values  = reshape(zeta, [1, 1, ctr.kmax]);
    zeta_diva = ones(ctr.imax, ctr.jmax, 1) .* k_values; % Use implicit expansion (broadcasting) to create the 3D matrix.

    % For some reason, there are weir infinities in the dudz, dvdz.
    % Those are in grounded ice and do not affect the simulation.
    % By definition: zeta(1)=0, zeta(ctr.kmax)=1;
    %zeta_m = flip(zeta);
    zeta_m = 1.0 - zeta_diva;

    % Eq. 36 Lipscomb et al. (2019).
    dudz = taux_diva .* H_diva .* zeta_m ./ eta_diva;
    dvdz = tauy_diva .* H_diva .* zeta_m ./ eta_diva;



    %for k = 1:ctr.kmax
        % In Kori, eta is defined as eta*H.
        %dudz(:,:,k) = taux .* H * ( 1.0 - zeta(k) ) ./ eta_diva(:,:,k);
        %dvdz(:,:,k) = tauy .* H * ( 1.0 - zeta(k) ) ./ eta_diva(:,:,k);
    %    dudz(:,:,k) = taux .* H * zeta_m(k) ./ eta_diva(:,:,k);
    %    dvdz(:,:,k) = tauy .* H * zeta_m(k) ./ eta_diva(:,:,k);
        
    %end


    %if any( isinf(dudz(:)) | isinf(dvdz(:)) )    
    %    error('Matrix contains Inf or -Inf values. Stopping execution.');
    %end


    % Strain rates including vertical derivatives.
    dudx_diva = repmat(dudx, [1, 1, ctr.kmax]);
    dvdy_diva = repmat(dvdy, [1, 1, ctr.kmax]);
    dudy_diva = repmat(dudy, [1, 1, ctr.kmax]);
    dvdx_diva = repmat(dvdx, [1, 1, ctr.kmax]);

    EffStr_diva = dudx_diva.^2 + dvdy_diva.^2 + ...
                    dudx_diva.*dvdy_diva + ...
                        0.25*(dudy_diva+dvdx_diva).^2 + ...
                            0.25 * (dudz.^2 + dvdz.^2);

    % Obtain 2D strain rates by vertically averaging the 3D DIVA strain rates.
    EffStr_diva = max(EffStr_diva,1e-12);
    EffStr      = mean(EffStr_diva, 3);

    %eta=0.5*H.*A.^(-1./par.n).*EffStr.^((1-par.n)/(2*par.n));
    %eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune

    % 3D strain rate.
    % Replicate ice rate factor A since it takes value from the base and not
    % vertical integration (Frank followed Ritz's approach).
    %H_diva = repmat(H, [1, 1, ctr.kmax]);
    A_diva = repmat(A, [1, 1, ctr.kmax]);

    %eta_diva=0.5*H.*A.^(-1./par.n).*EffStr_diva.^((1-par.n)/(2*par.n));
    eta_diva = 0.5 * H_diva .* A_diva.^(-1./par.n) ...
                    .* EffStr_diva.^((1-par.n)/(2*par.n));

    %eta_diva(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune

    % Vertical average of 3D eta_diva to pass it to 2D solver.
    eta = mean(eta_diva, 3);
    eta(MASK==0)=eta(MASK==0)./shelftune(MASK==0);  %VL: 2D shelftune


    %if any( isinf(eta_diva) )      
    %    error('Matrix contains Inf or -Inf values. Stopping execution.');
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    % Daniel: new grain size model.
    % Based on Austin and Evans (2007); Ranganathan et al. (2021).
    % Grain size model (built upon Ranganathan et al., 2021).
    dynamic_d_grain = true;
    if dynamic_d_grain == true
        d_grain = GrainSize(H,tmp,EffStr,zeta,ctr,par);

        %min(d_grain(:))
        %max(d_grain(:))
        
    % Constant grain size.
    else
        % Diffussion creep
        d_grain=5e-3; % tunable parameter? Bassis found low effect.
    end


    % Jablasco regularization implementation.
    % Bassis et al., (2021) regularization
    if ctr.bassis_reg==1
        % Glen flow law
        eta_glen=eta;
        eta_diff=H.*(A.^(-1/par.n))./(2*d_grain.^2);
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
        %eta(glMASK==6)=8.0e9; % Default: 1.0e7. Daniel: 1.0e10. Final: 8.0e9. Pseudo-transient: 0.5e5 or comment line.
        eta(glMASK==6)=mean(eta(glMASK==5)); % Take average of calving front viscosity
    end
    
end


