function [uxssa,uyssa,beta2,eta,dudx,dudy,dvdx,dvdy,su,ubx,uby,ux,uy, ...
    damage,NumStabVel,k,err]= ...
    SSAvelocity(ctr,par,su,Hmx,Hmy,gradmx,gradmy,signx,signy, ...
    uxssa,uyssa,H,HB,B,stdB,Asf,A,MASK,glMASK,HAF,HAFmx,HAFmy,cnt, ...
    nodeu,nodev,MASKmx,MASKmy,bMASK,uxsia,uysia,udx,udy,node,nodes, ...
    Mb,Melt,dtdx,dtdx2,VM,damage,ThinComp,shelftune)

% Kori-ULB
% Iterative solution to the SSA velocity (both pure SSA and hybrid model

    eps=1e-8;
    taudx=par.rho*par.g*Hmx.*sqrt(gradmx).*signx;
    taudy=par.rho*par.g*Hmy.*sqrt(gradmy).*signy;
    if ctr.uSSAexist==1 || cnt>1
        ussa=vec2h(uxssa,uyssa);    %VL: ussa on h-grid
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
    beta2=min(beta2,1e8);
    beta2(MASK==0)=0;
    betax=0.5*(beta2+circshift(beta2,[0 -1]));
    betay=0.5*(beta2+circshift(beta2,[-1 0]));

    if ctr.mismip>=1
        betax(:,1)=betax(:,2); % symmetric divide
        betax(1,:)=betax(3,:); % symmetry axis
        betax(ctr.imax,:)=betax(ctr.imax-2,:); % periodic BC
        betax(:,ctr.jmax)=0; % ocean
        betay(:,1)=betay(:,3); % symmetric divide
        betay(1,:)=betay(2,:); % symmetry axis
        betay(ctr.imax,:)=betay(ctr.imax-1,:); % periodic BC
        if ctr.mismip==2 % Thule setup
            betax(ctr.imax,:)=0;
            betay(ctr.imax,:)=0;
        end
    end

    % Viscosity on edge (or when viscosity cannot be calculated)
    if ctr.SSA>1 % hydrid model
        udx(HAFmx<0)=0; % no deformation for floating cells
        udy(HAFmy<0)=0;
    else % pure SSA
        udx=zeros(ctr.imax,ctr.jmax);
        udy=zeros(ctr.imax,ctr.jmax);
    end

    for ll=1:par.visciter % iteration over effective viscosity
        [eta,dudx,dvdy,dudy,dvdx]=EffVisc(A,uxssa,uyssa,H,par,MASK, ...
            glMASK,shelftune,ctr);

        % Jablasco damage.
        if ctr.damage==1 && cnt>1
            if ll==1
                if ctr.localdamage~=1
                    % compute transport from previous damage
                    if ctr.thinning==1 % Thinning component (very sensitive)
                        ThinComp = ThinningComponent(ctr,par,dudx,dvdy,dudy,dvdx,eta,H);
                        ThinComp(MASK==1)=0.0;
                        ThinComp(MASKlk==1)=0.0;
                    else
                        ThinComp = zeros(ctr.imax,ctr.jmax);
                    end
                    
                    dtr=TransportDamage(node,nodes,damage,Mb,Melt,ThinComp,H,glMASK,dtdx,dtdx2, ...
                                    uxssa,uyssa,ctr,cnt,bMASK,VM,par);
                else
                    dtr=zeros(ctr.imax,ctr.jmax);
                end
                
                % compute surface damage
                ds=SurfaceDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,MASK);
                % compute basal damage (and Kachuck term, necessary for transport)
                db=BasalDamageAlgorithms(ctr,par,dudx,dvdy,dudy,dvdx,eta,H,HAF);
                % Avoid basal damage on grounded ice
                db(MASK==1)=0.0;
                % Avoid damage on lakes
                db(MASKlk==1)=0.0;
                ds(MASKlk==1)=0.0;

                % total damage is sum of surface and basal damage
                % damage is limited to damlim
                dlim=max(0,min(db+ds,H.*par.dlim));
                damage=min(par.damlim*H,max(dlim,dtr));
                scale_eta=(H-min(damage,H-eps))./(H+eps);
            end
        else
            scale_eta=1;
            damage=zeros(ctr.imax,ctr.jmax);
        end
        eta=eta.*scale_eta;
        
        % PSEUDO-TRANSIENT METHOD.
        % First, implicit initialization to avoid zeros.
        if cnt < 0 % 0.1*ctr.nsteps, 20, 40, 200
            k=0.0;
            err=0.0;
            [uxs1,uys1,su,flagU,relresU,iterU]=SparseSolverSSA_daniel(nodeu,nodev, ...
                su,MASKmx,MASKmy,bMASK, ...
                H,eta,betax,betay,uxssa,uyssa,uxsia,uysia,udx,udy,taudx, ...
                taudy,ctr,par);

                uxssa=uxs1;
                uyssa=uys1;
                %---------iterative beta---------
                if cnt<=ctr.BetaIter
                    ussa=vec2h(uxssa,uyssa); %VL: ussa on h-grid
                    if ctr.u0>1e10
                        beta2=fg.*(ussa.^(1/ctr.m-1)).*Asf.^(-1/ctr.m);
                    else
                        beta2=fg.*(ussa.^(1/ctr.m-1)).*((ussa+ctr.u0).*Asf ...
                            /ctr.u0).^(-1/ctr.m);
                    end
                    beta2=min(beta2,1e8);
                    beta2(MASK==0)=0;
                    betax=0.5*(beta2+circshift(beta2,[0 -1]));
                    betay=0.5*(beta2+circshift(beta2,[-1 0]));
                    if ctr.mismip>=1
                        betax(:,1)=betax(:,2); % symmetric divide
                        betax(1,:)=betax(3,:); % symmetry axis
                        betax(ctr.imax,:)=betax(ctr.imax-2,:); % periodic BC
                        betax(:,ctr.jmax)=0; % ocean
                        betay(:,1)=betay(:,3); % symmetric divide
                        betay(1,:)=betay(2,:); % symmetry axis
                        betay(ctr.imax,:)=betay(ctr.imax-1,:); % periodic BC
                        if ctr.mismip==2 % Thule setup
                            betax(ctr.imax,:)=0;
                            betay(ctr.imax,:)=0;
                        end
                    end
                end
                %--------------------------------
        % Pseudo-transient solver.
        else
            [uxs1,uys1,k,err]=SolverSSA_pseudo_transient(nodeu,nodev, ...
                 su,MASKmx,MASKmy,bMASK,H,eta,uxssa,uyssa,...
                   uxsia,uysia,betax,betay,udx,udy,taudx,taudy,MASK,Asf,cnt,ctr,par);

        end


        % Pseudo-transient method.
        %[uxs1,uys1,betax,betay]=SparseSolverSSA_daniel_2(nodeu,nodev, ...
        %        su,MASKmx,MASKmy,bMASK,H,eta,uxssa,uyssa,...
        %            uxsia,uysia,udx,udy,taudx,taudy,MASK,Asf,cnt,ctr,par);
        

        duxs=sqrt((uxs1-uxssa).^2+(uys1-uyssa).^2);
        duxs(isnan(duxs))=0;
        uxssa=uxs1;
        uyssa=uys1;
        limit=sum(duxs(:))/(ctr.imax*ctr.jmax);

        if limit<par.visctol % Limit on convergence
            break;
        end
    end
    if ctr.NumCheck==1
        NumStabVel=[ll,limit,relresU,iterU,flagU];
    else
        NumStabVel=false;
    end
    uxssa(Hmx==0 & HAFmx>0)=0; % only for grounded ice sheet with H=0
    uyssa(Hmy==0 & HAFmy>0)=0;
    uxssa=min(max(-par.maxspeed,uxssa),par.maxspeed);
    uyssa=min(max(-par.maxspeed,uyssa),par.maxspeed);
    ux=uxssa;   %LZ2021
    uy=uyssa;   %LZ2021
    if ctr.SSAdiffus==2 % SSA as basal velocity only
        ux=ux+udx;
        uy=uy+udy;
    end
    if ctr.SSA==1 % basal sliding velocity SSA
        ubx=ux;
        uby=uy;
    else
        ubx=ux-udx;
        uby=uy-udy;
    end
end


