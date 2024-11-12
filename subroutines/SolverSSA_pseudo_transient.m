function [u,v,k,err]=SolverSSA_pseudo_transient(nodeu,nodev, ...
            s0,MASKmx,MASKmy,bMASK,H,eta,...
            u,v,betax,betay,usia,vsia,udx,udy,taudx,taudy,MASK,Asf,cnt,ctr,par)

% Kori-ULB
% Solving the SSA equation (both pure and hybrid SSA) applying the so-called
% "Pseudo-transient" method. This is a matrix-free approach that allows
% for extremely fast computation without strong memory requirements.
% Velocity solution for two two-dimensional ice shelf velocity field
% with kinematic boundary conditions.
% u-field on u-grid, viscosity on h-grid

        
    % Pseudo-time iteration parameters.
    rel   = 0.2;      % Relaxation between two consecutive solutions. 0.3
    %alpha = 1.0e-6;   % Pseudo-time step length. 1.0e-6. Critical to ensure convergence! Not increase too much.
    iter  = 200;       % Max number of iterations 500, 100. 50 also works.
    tol   = 2.0e-3;   % Tolerance to accept new solution. 1.0e-2, 1.0e-3 works. 1.0e-4 is better. 1.0e-2 is problematic.
    % Definitions for boundary conditions.
    A = 0.25*ctr.delta*par.rho*par.g*(1.-par.rho/par.rhow)*H.^2./eta;

    % Try good guess. Remember that in Kori eta = eta.*H.
    D = par.rho * H ./ ( 4.0 * eta * par.secperyear );
    eta_b = 0.5; % Ice bluk viscosity
    n_dim = 4.1; % Dimension factor in 2D.
    alpha = D .* (ctr.delta)^2 ./ ( (1.0+eta_b) * n_dim );

    
    % BOUNDARIES: j=1 and j=jmax.
    a = 2:ctr.imax;
    % BOUNDARIES: i=1 and i=imax.
    b = 2:ctr.jmax;


    dx_inv_2 = 1.0./(ctr.delta*ctr.delta);


    % Pseudo-time loop.
    err = 1.0;
    k = 0;
    while k < iter && err > tol

        u_old = u;
        v_old = v;

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
        if ctr.uSSAexist==1 || cnt>1
            ussa=vec2h(u,v);    %VL: ussa on h-grid
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
        %beta2=min(beta2,1e8);
        %beta2(MASK==0)=0;
        %betax=beta2;
        %betay=beta2;

        % Evaluate stress balance.
        stress_x = dx_inv_2 * ( fx_1 + fx_2 ) - betax .* u + taudx;
        stress_y = dx_inv_2 * ( fy_1 + fy_2 ) - betay .* v + taudy;

        u = u_old + alpha .* stress_x;    
        v = v_old + alpha .* stress_y;

        % BOUNDARY CONDITIONS IN VECTORIAL FORM.
        v_y = v - circshift(v,[1 0]); % (i-1,j)
        u_y = u - circshift(u,[1 0]); % (i-1,j)
        v_x = v - circshift(v,[0 1]); % (i,j-1)
        u_x = u - circshift(u,[0 1]); % (i,j-1)

        % BOUNDARIES: j=1 and j=jmax.
        % x-component.
        u(a,1)        = u(a,2) - 0.5 * ( - v_y(a,2) + A(a,2) );
        u(a,ctr.jmax) = u(a,ctr.jmax-1) + 0.5 * ( - v_y(a,ctr.jmax-1) + A(a,ctr.jmax-1) );

        % y-component.
        v(a,1)        = v(a,2) + u_y(a,2);
        v(a,ctr.jmax) = v(a,ctr.jmax-1) - u_y(a,ctr.jmax-1);
        

        % BOUNDARIES: i=1 and i=imax.
        % y-component.
        v(1,b)        = v(2,b) - 0.5 * ( - u_x(2,b) + A(2,b) );
        v(ctr.imax,b) = v(ctr.imax-1,b) + 0.5 * ( - u_x(ctr.imax-1,b) + A(ctr.imax-1,b) );

        % x-component.
        u(1,b)        = u(2,b) + v_x(2,b);
        u(ctr.imax,b) = u(ctr.imax-1,b) - v_x(ctr.imax-1,b);


        % Symmetry axis mismip quarter of a circle. Good!
        if ctr.mismip==2
            u(:,1) = - u(:,3);
            v(:,1) = v(:,3);

            u(1,:) = u(3,:);
            v(1,:) = -v(3,:);
        end

        % Not ready yet!
        if ctr.mismip==1
            u(:,1) = 0.0;
            %v(:,1) = 0.0;

            u(1,:) = u(3,:);
            %v(1,:) = -v(3,:); 

            v(ctr.imax,:) = 0.0;
            u(ctr.imax,:) = u(ctr.imax-1,:);
        end


        % At j=1:    du/dx = (-3*f[i+0]+4*f[i+1]-1*f[i+2])/(2*h) 
        % At j=jmax: du/dx = (1*f[i-2]-4*f[i-1]+3*f[i+0])/(2*h)

        %gamma_j1    = 0.5 * ( - v_y(a,2) + A(a,1) );
        %gamma_jmax = 0.5 * ( - v_y(a,ctr.jmax-1) + A(a,ctr.jmax) );

        %u(a,1)        = - ( 4.0 * u(a,2) - u(a,3) + 2.0 * gamma_j1 ) / 3.0;
        %u(a,ctr.jmax) = ( 4.0 * u(a,ctr.jmax-1) - u(a,ctr.jmax-2) + 2.0 * gamma_jmax ) / 3.0;

        %v(a,1)        = - ( 4.0 * v(a,2) - v(a,3) + 2.0 * u_y(a,2) ) / 3.0;
        %v(a,ctr.jmax) = ( 4.0 * v(a,ctr.jmax-1) - v(a,ctr.jmax-2) + 2.0 * u_y(a,ctr.jmax-1) ) / 3.0;

        % BOUNDARIES: i=1 and i=imax.
        %b = 2:ctr.jmax;
        
        % y-component.
        %gamma_i1    = 0.5 * ( - u_x(2,b) + A(1,b) );
        %gamma_imax = 0.5 * ( - u_x(ctr.imax-1,b) + A(ctr.imax-1,b) );

        %v(1,b)         = - ( 4.0 * v(2,b) - v(3,b) + 2.0 * gamma_i1 ) / 3.0;
        %v(ctr.imax,b) = ( 4.0 * v(ctr.imax-1,b) - v(ctr.imax-2,b) + 2.0 * gamma_imax ) / 3.0;

        %u(1,b)         = - ( 4.0 * u(2,b) - u(3,b) + 2.0 * v_x(2,b) ) / 3.0;
        %u(ctr.imax,b) = ( 4.0 * u(ctr.imax-1,b) - u(ctr.imax-2,b) + 2.0 * v_x(ctr.imax-1,b) ) / 3.0;


        % Update solution. Relaxation to avoid spurious results.
        u = rel * u_old + (1.0 - rel) * u;
        v = rel * v_old + (1.0 - rel) * v;

        % Update current error.
        dif = sqrt((u-u_old).^2 + (v-v_old).^2)./sqrt(u.^2 + v.^2);
        err = norm(dif,2);

        % Use difference for the step alpha.
        %dif_norm = dif / max(dif);
        %alpha = alpha_min * (1.0 - dif_norm) + alpha_max * dif_norm;

        k = k + 1;
    end


end


