function [u,v,betax,betay]=SparseSolverSSA_daniel_2(nodeu,nodev, ...
            s0,MASKmx,MASKmy,bMASK,H,eta,...
            u,v,usia,vsia,udx,udy,taudx,taudy,MASK,Asf,cnt,ctr,par)

% Kori-ULB
% Solving the SSA equation (both pure and hybrid SSA)
% sparse matrix solution for two two-dimensional ice shelf velocity field
% with kinematic boundary conditions
% u-field on u-grid, viscosity on h-grid
% u-velocities: quadrant U and Uv of solution matrix
% V-velocities: quadrant V and Vu of solution matrix
% Subsequent interleaving of u and v velocities (Quiquet et al., 2018) to
% improve stability and speed of the algorithm


    if cnt < 5

        ussa=vec2h(u,v); %VL: ussa on h-grid
        if par.ShelfPinning==1 && ctr.stdBexist==1 && ctr.inverse==0 
            fg=max(0,1-(HB-B)./stdB); % subgrid pinning points in ice shelf
        else
            fg=1;
        end
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

        limit=1.0e-5; % limit on effective viscosity gradients 1.0e-5. why??

        if ctr.SSAdiffus==2
            udx=zeros(ctr.imax,ctr.jmax);
            udy=zeros(ctr.imax,ctr.jmax);
        end

        R0=zeros(ctr.imax,ctr.jmax);
        U0=zeros(ctr.imax,ctr.jmax); % u(i,j)
        U1=zeros(ctr.imax,ctr.jmax); % u(i,j+1)
        U2=zeros(ctr.imax,ctr.jmax); % u(i,j-1)
        U3=zeros(ctr.imax,ctr.jmax); % u(i+1,j)
        U4=zeros(ctr.imax,ctr.jmax); % u(i-1,j)
        U5=zeros(ctr.imax,ctr.jmax); % u(i+1,j+1)
        U6=zeros(ctr.imax,ctr.jmax); % u(i+1,j-1)
        U7=zeros(ctr.imax,ctr.jmax); % u(i-1,j+1)
        U8=zeros(ctr.imax,ctr.jmax); % u(i-1,j-1)
        U9=zeros(ctr.imax,ctr.jmax); % periodic BC i=1
        U10=zeros(ctr.imax,ctr.jmax); % periodic BC i=imax
        Uv0=zeros(ctr.imax,ctr.jmax); % v(i,j)
        Uv1=zeros(ctr.imax,ctr.jmax); % v(i,j+1)
        Uv2=zeros(ctr.imax,ctr.jmax); % v(i-1,j)
        Uv3=zeros(ctr.imax,ctr.jmax); % v(i-1,j+1)

        % Frank.
        eta1=circshift(eta,[0 -1]); % eta(i,j+1)
        H1=circshift(H,[0 -1]); % H(i,j+1)
        
        % Daniel.
        %eta6=circshift(eta,[0 1]);
        %eta1=circshift(eta,[0 1]); % eta(i,j-1)
        %H1=circshift(H,[0 1]); % H(i,j+1)
        
        eta2=circshift(eta,[-1 0]); % eta(i+1,j)
        eta3=circshift(eta,[-1 -1]); % eta(i+1,j+1)
        eta4=circshift(eta,[1 0]); % eta(i-1,j)
        eta5=circshift(eta,[1 -1]); % eta(i-1,j+1)
        
        % Frank.
        dmudx=(eta1-eta)/ctr.delta;
        dmudy=0.25*(eta2+eta3-eta4-eta5)/ctr.delta;
        dmudx=min(limit,max(dmudx,-limit));
        dmudy=min(limit,max(dmudy,-limit));
        

        % Daniel.
        %dmudx=zeros(ctr.imax,ctr.jmax);
        %dmudy=zeros(ctr.imax,ctr.jmax);
        %dmudx=0.5.*(eta1-eta6)/ctr.delta;
        %dmudy=0.5*(eta+eta1-eta4-eta5)/ctr.delta;
        

        MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for SSA
        MASKb(2:ctr.imax-1,2:ctr.jmax-2)=1;
        MASKb(MASKmx==0)=0;

        %dmudx(MASKb==1)=0.0;

        % Frank.
        U0(MASKb==1)=-5.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
            betax(MASKb==1);
        U1(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)+ ...
            2.*dmudx(MASKb==1)/ctr.delta;
        U2(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
            2.*dmudx(MASKb==1)/ctr.delta;
        U3(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)+ ...
            0.5*dmudy(MASKb==1)/ctr.delta;
        U4(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
            0.5*dmudy(MASKb==1)/ctr.delta;
        Uv0(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
            (dmudx(MASKb==1)-0.5*dmudy(MASKb==1))/ctr.delta;
        Uv1(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
            (dmudx(MASKb==1)+0.5*dmudy(MASKb==1))/ctr.delta;
        Uv2(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
            (dmudx(MASKb==1)+0.5*dmudy(MASKb==1))/ctr.delta;
        Uv3(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
            (dmudx(MASKb==1)-0.5*dmudy(MASKb==1))/ctr.delta;


        % UNDER DEVELOPMENT.
        %U1(MASKb==1)=0.0;
        %U2(MASKb==1)=0.0;
        %U3(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2);
        %U4(MASKb==1)=0.0;

        % Daniel.
        % u(i,j)
        %U0(MASKb==1)=-5.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
        %    betax(MASKb==1);
        
        % u(i,j+1)
        %U1(MASKb==1)=4.*eta1(MASKb==1)/(ctr.delta^2);

        % u(i,j-1)
        %U2(MASKb==1)=4.*eta(MASKb==1)/(ctr.delta^2);

        % u(i+1,j)
        %U3(MASKb==1)=eta2(MASKb==1)/(ctr.delta^2);
        
        % u(i-1,j)
        %U4(MASKb==1)=eta(MASKb==1)/(ctr.delta^2);

        % v(i,j)
        %Uv0(MASKb==1)=-2.*eta(MASKb==1)/(ctr.delta^2) - eta2(MASKb==1)/(ctr.delta^2);
        
        % v(i,j+1)
        %Uv1(MASKb==1)=2.*eta1(MASKb==1)/(ctr.delta^2) + eta2(MASKb==1)/(ctr.delta^2);
        
        % v(i-1,j)
        %Uv2(MASKb==1)=3.*eta(MASKb==1)/(ctr.delta^2);
        
        % v(i-1,j+1)
        %Uv3(MASKb==1)=-2.*eta1(MASKb==1)/(ctr.delta^2) - eta(MASKb==1)/(ctr.delta^2);

        
        
        R0(MASKb==1)=-taudx(MASKb==1)-betax(MASKb==1).*udx(MASKb==1);

        MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for non-SSA
        MASKb(2:ctr.imax-1,2:ctr.jmax-2)=1;
        MASKb(MASKmx>0)=0;

        U0(MASKb==1)=1;
        R0(MASKb==1)=usia(MASKb==1);

        % boundary conditions

        if ctr.shelf==1 && ctr.mismip==0
            
            % j=1; contact with ocean (upwinding in x)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-1,1)=1;
            if ctr.basin==1
                U0(MASKb==1)=1;
                R0(MASKb==1)=u(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % FRANK.
            U0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta;
            U1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;
            Uv1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Uv3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
            
            % DANIEL.
            %U0(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
            %U1(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            %Uv1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            %Uv3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % DANIEL. Does not seem to produce changes if viscosity gradients have no limits.
            %U2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta; %-4.*eta(MASKb==1)/ctr.delta; % eta1 % u(i,j-1)
            %U0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            %Uv0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            %Uv2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % Frank.
            R0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.- ...
                par.rho/par.rhow);

            % Daniel.
            %R0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
            %    par.rho/par.rhow);

            % j=jmax-1; contact with ocean
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-1,ctr.jmax-1)=1; % ctr.jmax-1
            if ctr.basin==1
                U0(MASKb==1)=1;
                R0(MASKb==1)=u(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % FRANK.
            U0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            U2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
            Uv0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Uv2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % DANIEL.
            %U1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta; % u(i,j+1)
            %U0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta; % eta
            %Uv3(MASKb==1)=2.*eta1(MASKb==1)/ctr.delta; % Now: Uv1.With Uv3 gives nice results.
            %Uv2(MASKb==1)=-2.*eta1(MASKb==1)/ctr.delta;

            % DANIEL.
            %U1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta; % eta
            %U0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta; % eta
            %Uv1(MASKb==1)=2.*eta1(MASKb==1)/ctr.delta; 
            %Uv0(MASKb==1)=-2.*eta1(MASKb==1)/ctr.delta;


            % Frank.
            R0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
                par.rho/par.rhow);

            % Daniel.
            %R0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.- ...
            %    par.rho/par.rhow);

            
            % i=1; contact with ocean (v-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1,2:ctr.jmax-2)=1;
            if ctr.basin==1
                U0(MASKb==1)=1;
                R0(MASKb==1)=u(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % Frank.
            U0(MASKb==1)=-1/ctr.delta; % (i,j)
            U3(MASKb==1)=1./ctr.delta; % (i+1,j)
            Uv0(MASKb==1)=-1./ctr.delta; % (i,j)
            Uv1(MASKb==1)=1./ctr.delta; % (i,j+1)

            % Daniel. IT DOES NOT CHANGE ANYTHING,
            %U0(MASKb==1)=1/ctr.delta;
            %U4(MASKb==1)=-1./ctr.delta; % u(i-1,j)
            %Uv0(MASKb==1)=-1./ctr.delta;
            %Uv1(MASKb==1)=1./ctr.delta;

            R0(MASKb==1)=0;

            % i=imax; contact with ocean (v-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,2:ctr.jmax-2)=1;
            if ctr.basin==1
                U0(MASKb==1)=1;
                R0(MASKb==1)=u(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % FRANK.
            U0(MASKb==1)=1./ctr.delta;
            U4(MASKb==1)=-1./ctr.delta; % u(i-1,j)
            Uv2(MASKb==1)=-1./ctr.delta; % v(i-1,j)
            Uv3(MASKb==1)=1./ctr.delta; % v(i-1,j+1)

            % Daniel.
            %U3(MASKb==1)=1./ctr.delta;
            %U0(MASKb==1)=-1./ctr.delta;
            %Uv0(MASKb==1)=-1./ctr.delta;
            %Uv1(MASKb==1)=1./ctr.delta;

            R0(MASKb==1)=0;

            % j=jmax;
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(:,ctr.jmax)=1;
            U0(MASKb==1)=1;
            R0(MASKb==1)=NaN;

            % Model corners
            U0(1,1)=1;
            U1(1,1)=-1;
            U3(1,1)=-1;
            U5(1,1)=1;
            R0(1,1)=0;
            U0(1,ctr.jmax-1)=1;
            U2(1,ctr.jmax-1)=-1;
            U3(1,ctr.jmax-1)=-1;
            U6(1,ctr.jmax-1)=1;
            R0(1,ctr.jmax-1)=0;
            U0(ctr.imax,1)=1;
            U1(ctr.imax,1)=-1;
            U4(ctr.imax,1)=-1;
            U7(ctr.imax,1)=1;
            R0(ctr.imax,1)=0;
            U0(ctr.imax,ctr.jmax-1)=1;
            U2(ctr.imax,ctr.jmax-1)=-1;
            U4(ctr.imax,ctr.jmax-1)=-1;
            U8(ctr.imax,ctr.jmax-1)=1;
            R0(ctr.imax,ctr.jmax-1)=0;

        elseif ctr.shelf==1 && ctr.mismip>=1

            % j=1: ice divide (symmetric)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1:ctr.imax,1)=1;
            U0(MASKb==1)=1;
            U1(MASKb==1)=1;
            R0(MASKb==1)=0;

            % j=jmax-1; contact with ocean
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-1,ctr.jmax-1)=1;
            U0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            U2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
            Uv0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Uv2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
            R0(MASKb==1)=0.5*par.rho*par.g*(1.-par.rho/par.rhow)*H(MASKb==1).^2;

            % i=1: periodic boundary condition
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1,2:ctr.jmax)=1;
            U0(MASKb==1)=1;
            U9(MASKb==1)=-1;
            R0(MASKb==1)=0;

            if ctr.mismip==1
                % i=imax: periodic boundary condition
                MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
                MASKb(ctr.imax,2:ctr.jmax)=1;
                U0(MASKb==1)=1;
                U10(MASKb==1)=-1;
                R0(MASKb==1)=0;
            else
                % i=imax; contact with ocean (v-direction)
                MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
                MASKb(ctr.imax,2:ctr.jmax-2)=1;
                U0(MASKb==1)=1./ctr.delta;
                U4(MASKb==1)=-1./ctr.delta;
                Uv2(MASKb==1)=-1./ctr.delta;
                Uv3(MASKb==1)=1./ctr.delta;
                R0(MASKb==1)=0;
                U0(ctr.imax,ctr.jmax-1)=1;
                U2(ctr.imax,ctr.jmax-1)=-1;
                U4(ctr.imax,ctr.jmax-1)=-1;
                U8(ctr.imax,ctr.jmax-1)=1;
                R0(ctr.imax,ctr.jmax-1)=0;
            end

            % j=jmax;
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(:,ctr.jmax)=1;
            U0(MASKb==1)=1;
            R0(MASKb==1)=NaN;

        else
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1,1:ctr.jmax-1)=1;
            MASKb(ctr.imax,1:ctr.jmax-1)=1;
            MASKb(:,1)=1;
            MASKb(:,ctr.jmax-1)=1;
            U0(MASKb==1)=1;
            R0(MASKb==1)=0;

            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(:,ctr.jmax)=1;
            U0(MASKb==1)=1;
            R0(MASKb==1)=NaN;
        end

        % v-velocities: quadrant V and Vu of solution matrix

        S0=zeros(ctr.imax,ctr.jmax);
        V0=zeros(ctr.imax,ctr.jmax); % v(i,j)
        V1=zeros(ctr.imax,ctr.jmax); % v(i+1,j)
        V2=zeros(ctr.imax,ctr.jmax); % v(i-1,j)
        V3=zeros(ctr.imax,ctr.jmax); % v(i,j+1)
        V4=zeros(ctr.imax,ctr.jmax); % v(i,j-1)
        V5=zeros(ctr.imax,ctr.jmax); % v(i+1,j+1)
        V6=zeros(ctr.imax,ctr.jmax); % v(i-1,j+1)
        V7=zeros(ctr.imax,ctr.jmax); % v(i+1,j-1)
        V8=zeros(ctr.imax,ctr.jmax); % v(i-1,j-1)
        V9=zeros(ctr.imax,ctr.jmax); % Periodic BC on i=1
        V10=zeros(ctr.imax,ctr.jmax); % Periodic BC on i=imax-1
        V11=zeros(ctr.imax,ctr.jmax); % Symmetric ice divide v(i,j+3)
        Vu0=zeros(ctr.imax,ctr.jmax); % u(i,j)
        Vu1=zeros(ctr.imax,ctr.jmax); % u(i+1,j)
        Vu2=zeros(ctr.imax,ctr.jmax); % u(i,j-1)
        Vu3=zeros(ctr.imax,ctr.jmax); % u(i+1,j-1)

        % Frank.
        eta1=circshift(eta,[-1 0]); % eta(i+1,j)

        % Daniel
        %eta1=circshift(eta,[1 0]); % eta(i-1,j)
        
        % Frank.
        H1=circshift(H,[-1 0]); % H(i+1,j)

        % Daniel.
        %H1=circshift(H,[1 0]); % H(i-1,j)

        eta2=circshift(eta,[0 -1]); % eta(i,j+1)
        eta3=circshift(eta,[-1 -1]); % eta(i+1,j+1)
        % Frank.
        eta4=circshift(eta,[0 1]); % eta(i,j-1)
        
        eta5=circshift(eta,[-1 1]); % eta(i+1,j-1)
        
        % Frank.
        dmudy=(eta1-eta)/ctr.delta;
        dmudx=0.25*(eta2+eta3-eta4-eta5)/ctr.delta;
        dmudx=min(limit,max(dmudx,-limit));
        dmudy=min(limit,max(dmudy,-limit));

        % Daniel.
        %dmudy=(eta1-eta)/ctr.delta;
        %dmudx=(eta2-eta)/ctr.delta;


        MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for SSA
        MASKb(2:ctr.imax-2,2:ctr.jmax-1)=1; % 2:ctr.imax-2
        MASKb(MASKmy==0)=0;

        V0(MASKb==1)=-5.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2)- ...
            betay(MASKb==1); % -5.*
        V1(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)+ ...
            2.*dmudy(MASKb==1)/ctr.delta;
        V2(MASKb==1)=2.*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)- ...
            2.*dmudy(MASKb==1)/ctr.delta;
        V3(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)+ ...
            0.5*dmudx(MASKb==1)/ctr.delta;
        V4(MASKb==1)=0.5*(eta1(MASKb==1)+eta(MASKb==1))/(ctr.delta^2.)- ...
            0.5*dmudx(MASKb==1)/ctr.delta;
        Vu0(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
            (dmudy(MASKb==1)-0.5*dmudx(MASKb==1))/ctr.delta;
        Vu1(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)+ ...
            (dmudy(MASKb==1)+0.5*dmudx(MASKb==1))/ctr.delta;
        Vu2(MASKb==1)=1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
            (dmudy(MASKb==1)+0.5*dmudx(MASKb==1))/ctr.delta;
        Vu3(MASKb==1)=-1.5*(eta(MASKb==1)+eta1(MASKb==1))/(ctr.delta^2)- ...
            (dmudy(MASKb==1)-0.5*dmudx(MASKb==1))/ctr.delta;
        S0(MASKb==1)=-taudy(MASKb==1)-betay(MASKb==1).*udy(MASKb==1);

        MASKb=zeros(ctr.imax,ctr.jmax); % domain mask for non-SSA
        MASKb(2:ctr.imax-2,2:ctr.jmax-1)=1;
        MASKb(MASKmy>0)=0;
        V0(MASKb==1)=1;
        S0(MASKb==1)=vsia(MASKb==1);

        % boundary conditions

        if ctr.shelf==1 && ctr.mismip==0
            
            % i=1; contact with ocean
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1,2:ctr.jmax-1)=1;
            if ctr.basin==1
                V0(MASKb==1)=1;
                S0(MASKb==1)=v(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % FRANK.
            V0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta;
            V1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;
            Vu1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Vu3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % DANIEL.
            %eta1=circshift(eta,[1 0]); % [1 0] --> eta(i-1,j)
            %V2(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta; % v(i-1,j)
            %V0(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;  % v(i,j) 
            %Vu0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            %Vu2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % DANIEL 2. 
            %V0(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta; %1.*
            %V1(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            %Vu1(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            %Vu3(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % Frank.
            S0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.-par.rho/par.rhow);

            % Daniel.
            %S0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.-par.rho/par.rhow);

            % i=imax-1; contact with ocean
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax-1,2:ctr.jmax-1)=1; % ctr.imax-1
            if ctr.basin==1
                V0(MASKb==1)=1;
                S0(MASKb==1)=v(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            
            % FRANK.
            V0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
            V2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
            Vu0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
            Vu2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;

            % DANIEL.
            % Daniel: eta 4.
            %eta4=circshift(eta,[1 0]); % eta(i-1,j)
            %V1(MASKb==1)=4.*eta1(MASKb==1)/ctr.delta;
            %V0(MASKb==1)=-4.*eta1(MASKb==1)/ctr.delta;
            %Vu1(MASKb==1)=2.*eta1(MASKb==1)/ctr.delta;
            %Vu3(MASKb==1)=-2.*eta1(MASKb==1)/ctr.delta;

            % DANIEL 2.
            %V0(MASKb==1)=1.*eta1(MASKb==1)/ctr.delta;
            %V2(MASKb==1)=-1.*eta1(MASKb==1)/ctr.delta; % 4.0
            %Vu0(MASKb==1)=2.*eta1(MASKb==1)/ctr.delta;
            %Vu2(MASKb==1)=-2.*eta1(MASKb==1)/ctr.delta;

            % Frank.
            S0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.-par.rho/par.rhow);

            % Daniel.
            %S0(MASKb==1)=0.5*par.rho*par.g*H1(MASKb==1).^2.*(1.-par.rho/par.rhow);

            % j=1; contact with ocean (u-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-2,1)=1; % 2:ctr.imax-2,1
            if ctr.basin==1
                V0(MASKb==1)=1;
                S0(MASKb==1)=v(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            V0(MASKb==1)=-1/ctr.delta;
            V3(MASKb==1)=1./ctr.delta;
            Vu0(MASKb==1)=-1./ctr.delta;
            Vu1(MASKb==1)=1./ctr.delta;
            S0(MASKb==1)=0;


            % j=jmax; contact with ocean (u-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-2,ctr.jmax)=1; % 2:ctr.imax-2
            if ctr.basin==1
                V0(MASKb==1)=1;
                S0(MASKb==1)=v(MASKb==1);
                MASKb(bMASK==1)=0;
            end
            V0(MASKb==1)=1./ctr.delta;
            V4(MASKb==1)=-1./ctr.delta;
            Vu2(MASKb==1)=-1./ctr.delta;
            Vu3(MASKb==1)=1./ctr.delta;
            S0(MASKb==1)=0;

            % Model corners
            V0(1,1)=1;
            V1(1,1)=-1;
            V3(1,1)=-1;
            V5(1,1)=1;
            S0(1,1)=0;
            V0(ctr.imax-1,1)=1;
            V2(ctr.imax-1,1)=-1;
            V3(ctr.imax-1,1)=-1;
            V6(ctr.imax-1,1)=1;
            S0(ctr.imax-1,1)=0;
            V0(1,ctr.jmax)=1;
            V1(1,ctr.jmax)=-1;
            V4(1,ctr.jmax)=-1;
            V7(1,ctr.jmax)=1;
            S0(1,ctr.jmax)=0;
            V0(ctr.imax-1,ctr.jmax)=1;
            V2(ctr.imax-1,ctr.jmax)=-1;
            V4(ctr.imax-1,ctr.jmax)=-1;
            V8(ctr.imax-1,ctr.jmax)=1;
            S0(ctr.imax-1,ctr.jmax)=0;

            % i=imax;
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,:)=1;
            V0(MASKb==1)=1;
            S0(MASKb==1)=NaN;

        elseif ctr.shelf==1 && ctr.mismip>=1

            % j=1: ice divide (symmetric)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1:ctr.imax-1,1)=1;
            V0(MASKb==1)=1;
            V11(MASKb==1)=-1;
            S0(MASKb==1)=0;

            % j=jmax; contact with ocean (u-direction)
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(2:ctr.imax-2,ctr.jmax)=1;
            V0(MASKb==1)=1./ctr.delta;
            V4(MASKb==1)=-1./ctr.delta;
            Vu2(MASKb==1)=-1./ctr.delta;
            Vu3(MASKb==1)=1./ctr.delta;
            S0(MASKb==1)=0;

            % i=1: periodic boundary condition
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1,2:ctr.jmax)=1;
            V0(MASKb==1)=1;
            V9(MASKb==1)=1;
            S0(MASKb==1)=0;

            if ctr.mismip==1
                % i=imax-1: periodic boundary condition
                MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
                MASKb(ctr.imax-1,2:ctr.jmax)=1;
                V0(MASKb==1)=1;
                V10(MASKb==1)=1;
                S0(MASKb==1)=0;
            else
                % i=imax-1; contact with ocean
                MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
                MASKb(ctr.imax-1,2:ctr.jmax-1)=1;
                V0(MASKb==1)=4.*eta(MASKb==1)/ctr.delta;
                V2(MASKb==1)=-4.*eta(MASKb==1)/ctr.delta;
                Vu0(MASKb==1)=2.*eta(MASKb==1)/ctr.delta;
                Vu2(MASKb==1)=-2.*eta(MASKb==1)/ctr.delta;
                S0(MASKb==1)=0.5*par.rho*par.g*H(MASKb==1).^2.*(1.- ...
                    par.rho/par.rhow);
                V0(ctr.imax-1,ctr.jmax)=1;
                V2(ctr.imax-1,ctr.jmax)=-1;
                V4(ctr.imax-1,ctr.jmax)=-1;
                V8(ctr.imax-1,ctr.jmax)=1;
                S0(ctr.imax-1,ctr.jmax)=0;
            end

            % i=imax;
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,:)=1;
            V0(MASKb==1)=1;
            S0(MASKb==1)=NaN;

        else
            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(1:ctr.imax-1,1)=1;
            MASKb(1:ctr.imax-1,ctr.jmax)=1;
            MASKb(1,:)=1;
            MASKb(ctr.imax-1,:)=1;
            V0(MASKb==1)=1;
            S0(MASKb==1)=0;

            MASKb=zeros(ctr.imax,ctr.jmax); % boundary mask
            MASKb(ctr.imax,:)=1;
            V0(MASKb==1)=1;
            S0(MASKb==1)=NaN;
        end

        nodes=ctr.imax*ctr.jmax;
        V=[reshape(U0,nodes,1)
            reshape(V0,nodes,1)
            U1(U1~=0)
            U2(U2~=0)
            U3(U3~=0)
            U4(U4~=0)
            U5(U5~=0)
            U6(U6~=0)
            U7(U7~=0)
            U8(U8~=0)
            U9(U9~=0)
            U10(U10~=0)
            Uv0(Uv0~=0)
            Uv1(Uv1~=0)
            Uv2(Uv2~=0)
            Uv3(Uv3~=0)
            V1(V1~=0)
            V2(V2~=0)
            V3(V3~=0)
            V4(V4~=0)
            V5(V5~=0)
            V6(V6~=0)
            V7(V7~=0)
            V8(V8~=0)
            V9(V9~=0)
            V10(V10~=0)
            V11(V11~=0)
            Vu0(Vu0~=0)
            Vu1(Vu1~=0)
            Vu2(Vu2~=0)
            Vu3(Vu3~=0)
        ];

        row=[reshape(nodeu,nodes,1)
            reshape(nodev,nodes,1)
            nodeu(U1~=0)
            nodeu(U2~=0)
            nodeu(U3~=0)
            nodeu(U4~=0)
            nodeu(U5~=0)
            nodeu(U6~=0)
            nodeu(U7~=0)
            nodeu(U8~=0)
            nodeu(U9~=0)
            nodeu(U10~=0)
            nodeu(Uv0~=0)
            nodeu(Uv1~=0)
            nodeu(Uv2~=0)
            nodeu(Uv3~=0)
            nodev(V1~=0)
            nodev(V2~=0)
            nodev(V3~=0)
            nodev(V4~=0)
            nodev(V5~=0)
            nodev(V6~=0)
            nodev(V7~=0)
            nodev(V8~=0)
            nodev(V9~=0)
            nodev(V10~=0)
            nodev(V11~=0)
            nodev(Vu0~=0)
            nodev(Vu1~=0)
            nodev(Vu2~=0)
            nodev(Vu3~=0)
        ];
        nodeU1=circshift(nodeu,[0 -1]); %i,j+1
        nodeU2=circshift(nodeu,[0 1]); %i,j-1
        nodeU3=circshift(nodeu,[-1 0]); %i+1,j
        nodeU4=circshift(nodeu,[1 0]); %i-1,j
        nodeU5=circshift(nodeu,[-1 -1]); %i+1,j+1
        nodeU6=circshift(nodeu,[-1 1]); %i+1,j-1
        nodeU7=circshift(nodeu,[1 -1]); %i-1,j+1
        nodeU8=circshift(nodeu,[1 1]); %i-1,j-1
        nodeU9=circshift(nodeu,[-2 0]); % periodic BC i=1
        nodeU10=circshift(nodeu,[2 0]); %periodic BC i=imax
        nodeUv1=circshift(nodev,[0 -1]); %i,j+1
        nodeUv2=circshift(nodev,[1 0]); %i-1,j
        nodeUv3=circshift(nodev,[1 -1]); %i-1,j+1

        nodeV1=circshift(nodev,[-1 0]); %i+1,j
        nodeV2=circshift(nodev,[1 0]); %i-1,j
        nodeV3=circshift(nodev,[0 -1]); %i,j+1
        nodeV4=circshift(nodev,[0 1]); %i,j-1
        nodeV5=circshift(nodev,[-1 -1]); %i+1,j+1
        nodeV6=circshift(nodev,[1 -1]); %i-1,j+1
        nodeV7=circshift(nodev,[-1 1]); %i+1,j-1
        nodeV8=circshift(nodev,[1 1]); %i-1,j-1
        nodeV9=circshift(nodev,[-1 0]); % periodic BC i=1
        nodeV10=circshift(nodev,[1 0]); %periodic BC i=imax
        nodeV11=circshift(nodev,[0 -2]); % ice divide
        nodeVu1=circshift(nodeu,[-1 0]); %i+1,j
        nodeVu2=circshift(nodeu,[0 1]); %i,j-1
        nodeVu3=circshift(nodeu,[-1 1]); %i+1,j-1

        col=[reshape(nodeu,nodes,1)
            reshape(nodev,nodes,1)
            nodeU1(U1~=0)
            nodeU2(U2~=0)
            nodeU3(U3~=0)
            nodeU4(U4~=0)
            nodeU5(U5~=0)
            nodeU6(U6~=0)
            nodeU7(U7~=0)
            nodeU8(U8~=0)
            nodeU9(U9~=0)
            nodeU10(U10~=0)
            nodev(Uv0~=0)
            nodeUv1(Uv1~=0)
            nodeUv2(Uv2~=0)
            nodeUv3(Uv3~=0)
            nodeV1(V1~=0)
            nodeV2(V2~=0)
            nodeV3(V3~=0)
            nodeV4(V4~=0)
            nodeV5(V5~=0)
            nodeV6(V6~=0)
            nodeV7(V7~=0)
            nodeV8(V8~=0)
            nodeV9(V9~=0)
            nodeV10(V10~=0)
            nodeV11(V11~=0)
            nodeu(Vu0~=0)
            nodeVu1(Vu1~=0)
            nodeVu2(Vu2~=0)
            nodeVu3(Vu3~=0)
        ];

        % R=[reshape(R0,nodes,1)
        %     reshape(S0,nodes,1)];
        R(nodeu)=R0;
        R(nodev)=S0;
        R=R';
        R(isnan(R))=0;

        % construct sparse matrix
        A=sparse(row,col,V);

        % Cholesky factor and solve
        if ctr.ItSolv==1
            D=diag(diag(A));
            C1=tril(A);
            C2=D\triu(A);
            [s,flag,relres,iter]=bicgstab(A,R,par.veltol,par.veliter,C1,C2,s0);
            if flag>0
                s=A\R;
            end
        else
            s=A\R;
            [flag,relres,iter]=deal(false);
        end

        % Frank.
        u=s(nodeu);
        v=s(nodev);
    end



    if cnt >= 5
        % Pseudo-time iteration.
        rel = 0.25;
        alpha = 1.0e-6; % 1.0e-6. Critical to ensure convergence! Not increase too much.

        %fx_1 = zeros(ctr.imax,ctr.jmax);
        %fy_1 = zeros(ctr.imax,ctr.jmax);
        %fx_2 = zeros(ctr.imax,ctr.jmax);
        %fy_2 = zeros(ctr.imax,ctr.jmax);

        %fx = zeros(ctr.imax,ctr.jmax);
        %fy = zeros(ctr.imax,ctr.jmax);

        %M = zeros(ctr.imax,ctr.jmax);
        %M(4:ctr.imax-3,4:ctr.jmax-3) = 1;

        % Taudxy are in H grid.
        taudx = 0.5 * ( taudx + circshift(taudx, [0 -1]) );
        taudy = 0.5 * ( taudy + circshift(taudy, [-1 0]) );

        dx_inv_2 = 1.0./(ctr.delta*ctr.delta);

        iter = 100; % 50. Try 100
        err = 1.0;
        tol = 1.0e-3; % 1.0e-3 works. 1.0e-4 is better. 1.0e-2 can be problematic. 1.0e-6
        k = 0;

        while k < iter && err > tol

            u_old = u;
            v_old = v;

            % Loop over all points. Avoid last two due to boundary conditions.
%             for i=2:ctr.imax-1 % -1, -2
%                 for j=2:ctr.jmax-1 % -1, -2
%     
%                     % Handy definitions.
%                     u1 = u(i+1,j) + u(i+1,j-1);
%                     u2 = u(i,j) + u(i,j-1);
%                     u3 = u(i-1,j) + u(i-1,j-1);
%                     u4 = u(i,j+1) + u(i,j);
%                     u5 = u(i-1,j+1) + u(i-1,j);
% 
%                     v1 = v(i+1,j) + v(i,j);
%                     v2 = v(i+1,j-1) + v(i,j-1);
%                     v3 = v(i,j) + v(i-1,j);
%                     v4 = v(i,j-1) + v(i-1,j-1);
%                     v5 = v(i,j+1) + v(i-1,j+1);
% 
%                     phi = eta(i,j) * ( u2 - u3 + v3 - v4 );
% 
%                     % x-direction. 0.5 factor from averaging over two grid cells.
%                     fx_1(i,j) = 2.0 * ( eta(i,j+1) * ( 2.0 * ( u(i,j+1) - u(i,j) ) + v(i,j+1) - v(i-1,j+1) ) ...
%                                     - eta(i,j) * ( 2.0 * ( u(i,j) - u(i,j-1) ) + v(i,j) - v(i-1,j) ) );
% 
%                     fx_2(i,j) = 0.5 * ( eta(i+1,j) * ( u1 - u2 + v1 - v2 ) - phi );
% 
%                     % y-direction.
%                     fy_1(i,j) = 2.0 * ( eta(i+1,j) * ( 2.0 * ( v(i+1,j) - v(i,j) ) + u(i+1,j) - u(i+1,j-1) ) ...
%                                         - eta(i,j) * ( 2.0 * ( v(i,j) - v(i-1,j) ) + u(i,j) - u(i,j-1) ) );
% 
%                     fy_2(i,j) = 0.5 * ( eta(i,j+1) * ( v5 - v3 + u4 - u5 ) - phi );
% 
%                 end
%             end

            % Vectorial form.
            eta1=circshift(eta,[-1 0]); % (i+1,j)
            eta2=circshift(eta,[0 -1]); % (i,j+1)
    
            u1=circshift(u,[-1 0]); % (i+1,j)
            u2=circshift(u,[0 -1]); % (i,j+1)
            u3=circshift(u,[0 1]); % (i,j-1)
            u4=circshift(u,[1 0]); % (i-1,j)
            u5=circshift(u,[-1 1]); % (i+1,j-1)
            u6=circshift(u,[1 -1]); % (i-1,j+1)
            u7=circshift(u,[1 1]); % (i-1,j-1)
    
            v1=circshift(v,[-1 0]); % (i+1,j)
            v2=circshift(v,[0 -1]); % (i,j+1)
            v3=circshift(v,[0 1]); % (i,j-1)
            v4=circshift(v,[1 0]); % (i-1,j)
            v5=circshift(v,[-1 1]); % (i+1,j-1)
            v6=circshift(v,[1 -1]); % (i-1,j+1)
            v7=circshift(v,[1 1]); % (i-1,j-1)
    
            du1 = u1 + u5; % u(i+1,j) + u(i+1,j-1);
            du2 = u + u3;  % u(i,j) + u(i,j-1);
            du3 = u4 + u7; % u(i-1,j) + u(i-1,j-1);
            du4 = u2 + u;  % u(i,j+1) + u(i,j);
            du5 = u6 + u4; % u(i-1,j+1) + u(i-1,j);
    
            dv1 = v1 + v; % v(i+1,j) + v(i,j);
            dv2 = v5 + v3;  % v(i+1,j-1) + v(i,j-1);
            dv3 = v + v4; % v(i,j) + v(i-1,j);
            dv4 = v3 + v7;  % v(i,j-1) + v(i-1,j-1);
            dv5 = v2 + v6; % v(i,j+1) + v(i-1,j+1);
    
            phi = eta .* ( du2 - du3 + dv3 - dv4 );
            d3 = u - u3;
            d4 = v - v4;
    
            % x-direction.
            fx_1 = 2.0 * ( eta2 .* ( 2.0 * ( u2 - u ) + v2 - v6 ) ...
                          - eta .* ( 2.0 * d3 + d4 ) );
    
            fx_2 = 0.5 * ( eta1 .* ( du1 - du2 + dv1 - dv2 ) - phi );
    
            % y-direction.
            fy_1 = 2.0 * ( eta1 .* ( 2.0 * ( v1 - v ) + u1 - u5 ) ...
                          - eta .* ( 2.0 * d4 + d3 ) );
    
            fy_2 = 0.5 * ( eta2 .* ( dv5 - dv3 + du4 - du5 ) - phi );


            % We need to calculate beta herein to update it with the new velocity field.
            ussa=vec2h(u,v);
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
            % Frank.
            betax=0.5*(beta2+circshift(beta2,[0 -1]));
            betay=0.5*(beta2+circshift(beta2,[-1 0]));


            % Without stagerring results are good already! Run long sims to compare.
            % Perhaps beta should not be stagerred to avoid higher vel on second and thirth quadrants.
            % Stager fx and fy. TRY NOW THAT IT Is STABLE! 

            % Issues with BC
            %fx_1 = 0.5 * ( fx_1 + circshift(fx_1, [0 -1]) );
            %fx_2 = 0.5 * ( fx_2 + circshift(fx_2, [0 -1]) );
            %fy_1 = 0.5 * ( fy_1 + circshift(fy_1, [-1 0]) );
            %fy_2 = 0.5 * ( fy_2 + circshift(fy_2, [-1 0]) );
            
            % No issues with BC but same higher vel on second and thirth quadrants.
            %fx_1 = 0.5 * ( fx_1 + circshift(fx_1, [0 1]) );
            %fx_2 = 0.5 * ( fx_2 + circshift(fx_2, [0 1]) );
            %fy_1 = 0.5 * ( fy_1 + circshift(fy_1, [1 0]) );
            %fy_2 = 0.5 * ( fy_2 + circshift(fy_2, [1 0]) );

            % Evaluate stress balance.
            fx = dx_inv_2 * ( fx_1 + fx_2 ) - betax .* u;
            fy = dx_inv_2 * ( fy_1 + fy_2 ) - betay .* v;


            %fx_1 = ( fx_1 + circshift(fx_1, [0 -1]) + circshift(fx_1, [0 1]) ) / 3.0;
            %fx_2 = ( fx_2 + circshift(fx_2, [0 -1]) + circshift(fx_2, [0 1]) ) / 3.0;
            %fy_1 = ( fy_1 + circshift(fy_1, [-1 0]) + circshift(fy_1, [1 0]) ) / 3.0;
            %fy_2 = ( fy_2 + circshift(fy_2, [-1 0]) + circshift(fy_2, [1 0]) ) / 3.0;
            
            %fx(MASK==1) = dx_inv_2 * ( fx_1(MASK==1) + fx_2(MASK==1) ) - betax(MASK==1) .* u(MASK==1);
            %fy(MASK==1) = dx_inv_2 * ( fy_1(MASK==1) + fy_2(MASK==1) ) - betay(MASK==1) .* v(MASK==1);
            
            %fx(M==1) = dx_inv_2 * ( fx_1(M==1) + fx_2(M==1) ) - betax(M==1) .* u(M==1);
            %fy(M==1) = dx_inv_2 * ( fy_1(M==1) + fy_2(M==1) ) - betay(M==1) .* v(M==1);

            % New velocity solution. Correct sing?
            % Sign convention: SSA stress balance in Frank's.
            % R0(MASKb==1)=-taudx(MASKb==1)-betax(MASKb==1).*udx(MASKb==1);
            u = u_old + alpha * ( fx - taudx );
            v = v_old + alpha * ( fy - taudy );
            

            % Definitions for Boundary conditions.
            A = 0.25*ctr.delta*par.rho*par.g*(1.-par.rho/par.rhow)*H.^2./eta;

            % BOUNDARY CONDITIONS IN VECTORIAL FORM.
            v_y = v - circshift(v,[1 0]); % (i-1,j)
            u_y = u - circshift(u,[1 0]); % (i-1,j)
            v_x = v - circshift(v,[0 1]); % (i,j-1)
            u_x = u - circshift(u,[0 1]); % (i,j-1)
            
            % Boundaries: j=1 and j=jmax.
            a = 2:ctr.imax;
            
            % x-component.
            u(a,1)        = u(a,2) - 0.5 * ( - v_y(a,2) + A(a,2) );
            u(a,ctr.jmax) = u(a,ctr.jmax-1) + 0.5 * ( - v_y(a,ctr.jmax-1) + A(a,ctr.jmax-1) );
    
            % y-component.
            v(a,1)        = v(a,2) + u_y(a,2);
            v(a,ctr.jmax) = v(a,ctr.jmax-1) - u_y(a,ctr.jmax-1);
            
    
            % Boundaries: i=1 and i=imax.
            b = 2:ctr.jmax;
            
            % y-component.
            v(1,b)        = v(2,b) - 0.5 * ( - u_x(2,b) + A(2,b) );
            v(ctr.imax,b) = v(ctr.imax-1,b) + 0.5 * ( - u_x(ctr.imax-1,b) + A(ctr.imax-1,b) );
    
            % x-component.
            u(1,b)        = u(2,b) + v_x(2,b);
            u(ctr.imax,b) = u(ctr.imax-1,b) - v_x(ctr.imax-1,b);



            % Boundaries: j=1 and j=jmax.
%             for i=2:ctr.imax
%             %for i=2:ctr.imax-1
%             
%                 % x-component. Should it be A(i,1) instead and analogously everywhere else?
%                 v_y = v(i,2) - v(i-1,2);
%                 u(i,1) = u(i,2) - 0.5 * ( - v_y + A(i,1) );
% 
%                 v_y = v(i,ctr.jmax-1) - v(i-1,ctr.jmax-1);
%                 u(i,ctr.jmax) = u(i,ctr.jmax-1) + 0.5 * ( - v_y + A(i,ctr.jmax-1) );
% 
%                 % y-component.
%                 v(i,1)        = v(i,2) + ( u(i,2) - u(i-1,2) );
%                 v(i,ctr.jmax) = v(i,ctr.jmax-1) - ( u(i,ctr.jmax-1) - u(i-1,ctr.jmax-1) );
% 
%             end
% 
%             % Boundaries: i=1 and i=imax.
%             for j=2:ctr.jmax
%             %for j=2:ctr.jmax-1
%             
%                 % y-component.
%                 u_x = u(2,j) - u(2,j-1);
%                 v(1,j) = v(2,j) - 0.5 * ( - u_x + A(1,j) );
%                 %v(1,j) = v(2,j) - 0.5 * ( - u_x + 0.5*(A(1,j)+A(2,j)) );
%                 
%                 u_x = u(ctr.imax-1,j) - u(ctr.imax-1,j-1);
%                 v(ctr.imax,j) = v(ctr.imax-1,j) + 0.5 * ( - u_x + A(ctr.imax-1,j) );
%                 %v(ctr.imax,j) = v(ctr.imax-1,j) + 0.5 * ( - u_x + 0.5*(A(ctr.imax-1,j)+A(ctr.imax-2,j)) );
% 
%                 % x-component.
%                 u(1,j)        = u(2,j) + ( v(2,j) - v(2,j-1) );
%                 u(ctr.imax,j) = u(ctr.imax-1,j) - ( v(ctr.imax-1,j) - v(ctr.imax-1,j-1) );
% 
%             end


            % Error.
            dif = sqrt((u-u_old).^2 + (v-v_old).^2)./sqrt(u.^2 + v.^2);
            %err = max(dif,[],'all');
            err = norm(dif,2);

            % Update solution. Relaxation to avoid spurious results.
            u = rel * u_old + (1.0 - rel) * u;
            v = rel * v_old + (1.0 - rel) * v;

            % Bounds?
            %if cnt < 200 % 300 for u<3.0
            %    u(u<3.0) = 3.0;  % 5.0
            %    v(v<3.0) = 3.0;  % 5.0
            %end

            u(u>5.0e3) = 5.0e3;
            v(v>5.0e3) = 5.0e3;

            k = k + 1;
        end
    end


end


