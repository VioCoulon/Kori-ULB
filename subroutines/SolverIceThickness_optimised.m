function H=SolverIceThickness_optimised(Mb,H,u,v,MASK,d,B,SLR,par,ctr)

% Kori-ULB
% Optimized explicit solver of the ice thickness equation

    % Daniel's explicit calculation of ice thickness. Factor 0.75.
    dt_CFL = 0.5 * ctr.delta ./ ( abs(u) + abs(v) );
    dt_CFL = min(dt_CFL);


    if ctr.dt > dt_CFL
        fprintf('\n Set timestep too short for optimization.');
        fprintf('\n Current dt   = %12.2f\n\n', ctr.dt);
        fprintf('\n dt_CFL limit = %12.2f\n\n', dt_CFL);
    end

    %dt = min(ctr.dt,dt_CFL);
    dtdx=ctr.dt/ctr.delta;

    %dtdx=ctr.dt/(ctr.delta*ctr.delta);
    %dtdx2=ctr.dt/(2.*ctr.delta);

    R0 = Mb*ctr.dt;
    % This cannot be according to CalvMIP Exp 1!
    %R0(MASK==0)=0.0;
    %R0=zeros(ctr.imax,ctr.jmax);

    % Staggered definitions.
    H1=circshift(H,[-1 0]); % (i+1,j)
    H2=circshift(H,[0 -1]); % (i,j+1)
    %H3=circshift(H,[0 1]); % (i,j-1)
    %H4=circshift(H,[1 0]); % (i-1,j)

    %H5=circshift(H,[-2 0]); % (i+2,j)
    %H6=circshift(H,[2 0]); % (i-2,j)
    %H7=circshift(H,[0 -2]); % (i,j+2)
    %H8=circshift(H,[0 2]); % (i,j-2)
    
    %H9=circshift(H,[-3 0]); % (i+3,j)
    %H10=circshift(H,[3 0]); % (i-3,j)
    %H11=circshift(H,[0 -3]); % (i,j+3)
    %H12=circshift(H,[0 3]); % (i,j-3)

    
    %q_x     = zeros(ctr.imax,ctr.jmax);
    %q_y     = zeros(ctr.imax,ctr.jmax);
    %qx_grad = zeros(ctr.imax,ctr.jmax);
    %qy_grad = zeros(ctr.imax,ctr.jmax);
    
    % Boolean indexes.
    %up = u>0.0;
    %um = u<0.0;
    %vp = v>0.0;
    %vm = v<0.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diffusivity correction?
    % d has velocity units.
    %d1 = circshift(d, [1 0]); % (i-1,j)
    %d2 = circshift(d, [0 1]); % (i,j-1)
    %d3 = circshift(d, [1 1]); % (i-1,j-1)
    %d_u = 0.5 * ( d + d1 );
    %d_v = 0.5 * ( d + d2 );
    %s_u = sign(u);
    %s_v = sign(v);

    %M = zeros(ctr.imax,ctr.jmax); %dMASK (floating=0, grounded=1)
    %M( MASK==1 | MASK==2 ) = 1;

    %M1 = M + ( 1 - M ) * ( 1 - par.rho / par.rhow );
    %M2 = B .* M + SLR .* ( 1 - M );

    %u = u - 2.0 * s_u .* d_u .* MASKfac1;
    %v = v - 2.0 * s_v .* d_v .* MASKfac1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSISTENT WITH FRANK.
    %d1=circshift(d,[0 1]); % d(i,j-1)
    %d2=circshift(d,[1 0]); % d(i-1,j)
    %d3=circshift(d,[1 1]); % d(i-1,j-1)

    %H1 = circshift(H,[0 1]); % d(i,j-1)
    %%H2 = circshift(H,[1 0]); % d(i-1,j)
    %H3 = circshift(H,[0 -1]); % d(i,j+1)
    %H4 = circshift(H,[-1 0]);  % d(i+1,j)

    %u1 = circshift(u,[0 1]); % d(i,j-1)
    %v2 = circshift(v,[1 0]); % d(i-1,j)

    % Diffusivity terms.
    %d_u  = 2.0 * ( d + d2 );
    %d_u1 = 2.0 * ( d1 + d3 );
    
    %d_v  = 2.0 * ( d + d1 );
    %d_v2 = 2.0 * ( d2 + d3 );

    % Terms on each direction.
    %f_x = ( u - d_u .* M1 ) .* H3 - ( u1 + d_u1 .* M1 ) .* H1;
    %f_y = ( v - d_v .* M1 ) .* H4 - ( v2 + d_v2 .* M1 ) .* H2;

    % Diffusivity staggered onto H-grid.
    %D = d + d1 + d2 + d3;
    
    % Centred term on H_ij.
    %F = ( u - u1 + v - v2 + 4.0 * D .* M1 ) .* H;

    % Inhomogeneous diff.
    %M2_1 = circshift(M2,[0 -1]);
    %M2_2 = circshift(M2,[0 1]);
    %M2_3 = circshift(M2,[-1 0]);
    %M2_4 = circshift(M2,[1 0]);

    % Diffusivity in inhomogeneous term.
    %D0 = - 4.0*M2.*D + M2_1.*d_u + M2_2.*d_u1 + M2_3.*d_v + M2_4.*d_v2 ;

    %T0 = - H .* V0 - H3 .* V1 - H1 .* V2 - H4 .* V3 - H2 .* V4;

    %H_new = H + 0.5 * dtdx * ( f_x + f_y + F + D0 + T0) + R0;
    %H_new = H + 0.5 * dtdx * ( f_x + f_y + F ) + R0;




%     d1=circshift(d,[0 1]); % d(i,j-1)
%     d2=circshift(d,[1 0]); % d(i-1,j)
%     d3=circshift(d,[1 1]); % d(i-1,j-1)
% 
%     dMASK=zeros(ctr.imax,ctr.jmax); %dMASK (floating=0, grounded=1)
%     dMASK(MASK==1 | MASK==2)=1;
% 
%     um1=circshift(u,[0 1]); % u(i,j-1)
%     vm1=circshift(v,[1 0]); % v(i-1,j)
% 
%     dipx=dtdx*(d+d2);
%     dimx=dtdx*(d1+d3);
%     dipy=dtdx*(d+d1);
%     dimy=dtdx*(d2+d3);
% 
%     MASKfac1=dMASK+(1-dMASK)*(1-par.rho/par.rhow);
%     MASKfac2=B.*dMASK+SLR.*(1-dMASK);
% 
%     V0=2*(d+d1+d2+d3).*MASKfac1*dtdx+dtdx2*(u-um1+v-vm1); % i,j
%     V1=-dipx.*circshift(MASKfac1,[0 -1])+dtdx2*u; % i,j+1
%     V2=-dimx.*circshift(MASKfac1,[0 1])-dtdx2*um1; % i,j-1
%     V3=-dipy.*circshift(MASKfac1,[-1 0])+dtdx2*v; % i+1,j
%     V4=-dimy.*circshift(MASKfac1,[1 0])-dtdx2*vm1; % i-1,j

    %V0(MASK==0)=0; % note that for shelf=1, MASK=glMASK in the call
    %V1(MASK==0)=0;
    %V2(MASK==0)=0;
    %V3(MASK==0)=0;
    %V4(MASK==0)=0;

%     MASKb=zeros(ctr.imax,ctr.jmax);
%     MASKb(1,:)=1;
%     MASKb(ctr.imax,:)=1;
%     MASKb(:,1)=1;
%     MASKb(:,ctr.jmax)=1;
%     V1(MASKb==1)=0;
%     V2(MASKb==1)=0;
%     V3(MASKb==1)=0;
%     V4(MASKb==1)=0;
% 
%     H1 = circshift(H,[0 -1]);
%     H2 = circshift(H,[0 1]);
%     H3 = circshift(H,[-1 0]);
%     H4 = circshift(H,[1 0]);
% 
%     M1 = circshift(MASKfac2,[0 -1]);
%     M2 = circshift(MASKfac2,[0 1]);
%     M3 = circshift(MASKfac2,[-1 0]);
%     M4 = circshift(MASKfac2,[1 0]);
% 
%     % Grounded.
%     H_new = Mb*ctr.dt + H - H.*V0 - H1.*V1 - H2.*V2 - H3.*V3 - H4.*V4 - ...
%             MASKfac2.*(d+d1+d2+d3)*2*dtdx + M1.*dipx + ...
%             M2.*dimx + M3.*dipy+ M4.*dimy;
% 
%     % Reversed sign?
%     H_new = Mb*ctr.dt + H - H.*V0 - H1.*V1 - H2.*V2 - H3.*V3 - H4.*V4 + ...
%                 MASKfac2.*(d+d1+d2+d3)*2*dtdx - M1.*dipx - ...
%                 M2.*dimx - M3.*dipy - M4.*dimy;

    % Floating regions.
    %H_new(MASK==0) = Mb(MASK==0)*ctr.dt + H(MASK==0) ...
    %                - H(MASK==0).*dtdx2.*(u(MASK==0)-um1(MASK==0)+v(MASK==0)-vm1(MASK==0)) ...
    %                - H1(MASK==0).*dtdx2.*u(MASK==0) ...
    %                + H2(MASK==0).*dtdx2.*um1(MASK==0) ...
    %                - H3(MASK==0).*dtdx2.*v(MASK==0) ...
    %                + H4(MASK==0).*dtdx2.*vm1(MASK==0) ...
    %                + M1(MASK==0).*dipx(MASK==0) + M2(MASK==0).*dimx(MASK==0) ...
    %                + M3(MASK==0).*dipy(MASK==0)+ M4(MASK==0).*dimy(MASK==0);

    %H_new(MASK==0) = Mb(MASK==0)*ctr.dt + H(MASK==0) - ...
    %                    MASKfac2(MASK==0).*(d(MASK==0)+d1(MASK==0)+d2(MASK==0)+d3(MASK==0))*2*dtdx ...
    %                    + M1(MASK==0).*dipx(MASK==0) + ...
    %                    M2(MASK==0).*dimx(MASK==0) ...
    %                    + M3(MASK==0).*dipy(MASK==0) ...
    %                    + M4(MASK==0).*dimy(MASK==0);
    
                        % Staggered velocities.
    %u1 = circshift(u, [0 1]); % (i,j-1)
    %v1 = circshift(v, [1 0]); % (i-1,j)

    % Ice fluxes.
    %q_x(up) = u(up) .* H(up);
    %q_x(um) = u1(um) .* H(um);

    %q_y(vp) = v(vp) .* H(vp);
    %q_y(vm) = v1(vm) .* H(vm);

    % Satggered fluxes.
    %q_x1 = circshift(q_x, [0 1]); % (i,j-1)
    %q_x2 = circshift(q_x, [0 -1]); % (i,j+1)

    %q_y1 = circshift(q_y, [1 0]); % (i-1,j)
    %q_y2 = circshift(q_y, [-1 0]); % (i+1,j)

    % Fluxed horizontal gradients.
    %qx_grad(up) = q_x(up) - q_x1(up);
    %qx_grad(um) = q_x2(um) - q_x(um);

    %qy_grad(vp) = q_y(vp) - q_y1(vp);
    %qy_grad(vm) = q_y2(vm) - q_y(vm);

    % New ice thickness.
    %H_new = H - dtdx * ( qx_grad + qy_grad ) + R0;

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set velocities at the bordes to zero?
    % Use fluxes instead.
    Hx_stag = 0.5 * ( H + H2 ); % (i,j+1)
    Hy_stag = 0.5 * ( H + H1 ); % (i+1,j)

    q_x = u .* Hx_stag;
    q_y = v .* Hy_stag;

    qx_grad = q_x - circshift(q_x, [0 1]); % (i,j-1)
    qy_grad = q_y - circshift(q_y, [1 0]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dffusive correction? IT GIVES ALL ZEROS FOR PURELY SSA.
    %d1=circshift(d,[0 1]); % d(i,j-1)
    %d2=circshift(d,[1 0]); % d(i-1,j)   
    %d3=circshift(d,[1 1]); % d(i-1,j-1)

    %dMASK=zeros(ctr.imax,ctr.jmax); %dMASK (floating=0, grounded=1)
    %dMASK(MASK==1 | MASK==2)=1;

    %dipx = 0.5 * dtdx*(d+d2);
    %dimx = 0.5 * dtdx*(d1+d3);
    %dipy = 0.5 * dtdx*(d+d1);
    %dimy = 0.5 * dtdx*(d2+d3);
    %MASKfac1=dMASK+(1-dMASK)*(1-par.rho/par.rhow);
    %MASKfac2 = B .* dMASK + SLR .* (1-dMASK);

    %D = - MASKfac2.*(d+d1+d2+d3) * dtdx + ...
    %    circshift(MASKfac2,[0 -1]) .* dipx + ...
    %    circshift(MASKfac2,[0 1]) .* dimx + ...
    %    circshift(MASKfac2,[-1 0]) .* dipy + ...
    %    circshift(MASKfac2,[1 0]) .* dimy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    % New ice thickness.
    H_new = H - dtdx * ( qx_grad + qy_grad ) + R0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %H_new(:,1)        = H(:,1) + dtdx * ...
    %                            ( q_x(:,3) - q_x(:,2) + qy_grad(:,2) ) + R0(:,1);
    %H_new(:,ctr.jmax) = H(:,ctr.jmax) + dtdx * ...
    %                            ( q_x(:,ctr.jmax-1) - q_x(:,ctr.jmax-2) + qy_grad(:,ctr.jmax-1) ) + R0(:,ctr.jmax);

    %H_new(1,:)        = H(1,:) + dtdx * ...
    %                            ( qx_grad(2,:) + q_y(3,:) - q_y(2,:) ) + R0(1,:);
    %H_new(ctr.imax,:) = H(ctr.imax,:) + dtdx * ...
    %                            ( qx_grad(ctr.imax-1,:) + q_y(ctr.imax-2,:) - q_y(ctr.imax-1,:) ) + R0(ctr.imax,:);




    % Calving front.
    %H1=circshift(H,[-1 0]); % (i+1,j)
    %H2=circshift(H,[0 -1]); % (i,j+1)
    %H3=circshift(H,[0 1]); % (i,j-1)
    %H4=circshift(H,[1 0]); % (i-1,j)


    % ALL CASES.
    %u1 = circshift(u_stag,[0 1]); % (i,j-1)
    %v1 = circshift(v_stag,[1 0]); % (i-1,j)

    %H_new(calv & up & vp) = H(calv & up & vp) + dtdx * ...
    %                            ( u1(calv & up & vp) .* ( H3(calv & up & vp) - H(calv & up & vp) ) +...
    %                              v1(calv & up & vp) .* ( H4(calv & up & vp)  - H(calv & up & vp) ) ) +...
    %                            R0(calv & up & vp);

    %H_new(calv & um & vp) = H(calv & um & vp) + dtdx * ...
    %                                ( u(calv & um & vp) .* ( H(calv & um & vp) - H2(calv & um & vp) ) +...
    %                                  v1(calv & um & vp) .* ( H4(calv & um & vp)  - H(calv & um & vp) ) ) +...
    %                                R0(calv & um & vp);
    
    %H_new(calv & up & vm) = H(calv & up & vm) + dtdx * ...
    %                                    ( u1(calv & up & vm) .* ( H3(calv & up & vm) - H(calv & up & vm) ) +...
    %                                      v(calv & up & vm) .* ( H(calv & up & vm)  - H1(calv & up & vm) ) ) +...
    %                                    R0(calv & up & vm);
                                                                             
    %H_new(calv & um & vm) = H(calv & um & vm) + dtdx * ...
    %                                        ( u(calv & um & vm) .* ( H(calv & um & vm) - H2(calv & um & vm) ) +...
    %                                          v(calv & um & vm) .* ( H(calv & um & vm)  - H1(calv & um & vm) ) ) +...
    %                                        R0(calv & um & vm);
    
    

    % Symmetric. INTERESTING
    %H_x(up) = - 0.5 * ( H2(up) - H3(up) );
    %H_x(um) = 0.5 * ( H3(um) - H2(um) );
    %H_y(vp) = - 0.5 * ( H1(vp) - H4(vp) );
    %H_y(vm) = 0.5 * ( H4(vm) - H1(vm) );

    % BEST SO FAR.
    %H_x(up) = - ( H2(up) - H(up) );
    %H_x(um) =  H3(um) - H(um);
    %H_y(vp) = - ( H1(vp) - H(vp) );
    %H_y(vm) = H4(vm) - H(vm);

    % Second order derivatives (upstream).
    %H_x(up) = - 0.5 * ( H(up) - 4.0 * H3(up) + 3.0 * H8(up));
    %H_x(um) = 0.5 * ( H(um) - 4.0 * H2(um) + 3.0 * H7(um) );
    %H_y(vp) = - 0.5 * ( H(vp) - 4.0 * H4(vp) + 3.0 * H6(vp));
    %H_y(vm) = 0.5 * ( H(vm) - 4.0 * H1(vm) + 3.0 * H5(vm) );

    % Second order derivatives (downstream).
    %H_x(up) = - 0.5 * ( H7(up) - 4.0 * H2(up) + 3.0 * H(up));
    %H_x(um) = 0.5 * ( H8(um) - 4.0 * H3(um) + 3.0 * H(um) );
    %H_y(vp) = - 0.5 * ( H5(vp) - 4.0 * H1(vp) + 3.0 * H(vp));
    %H_y(vm) = 0.5 * ( H6(vm) - 4.0 * H4(vm) + 3.0 * H(vm) );

    % Try third order. (1 point downstream)
    % f_x = (1*f[i-2]-6*f[i-1]+3*f[i+0]+2*f[i+1]) / 6.0;
    %a = 1.0 / 6.0;
    %H_x(up) = - a * ( 2.0 * H2(up) + 3.0 * H(up) - 6.0 * H3(up) + H2(up) );
    %H_x(um) = a * ( 2.0 * H3(um) + 3.0 * H(um) - 6.0 * H2(um) + H7(um) );
    %H_y(vp) = - a * ( 2.0 * H1(vp) + 3.0 * H(vp) - 6.0 * H4(vp) + H6(vp) );
    %H_y(vm) = a * ( 2.0 * H4(vm) + 3.0 * H(vm) - 6.0 * H1(vm) + H5(vm) );

    % Third order. (all points upstream)
    % (-2*f[i-3]+9*f[i-2]-18*f[i-1]+11*f[i+0])/(6*1.0*h**1)
    %a = 1.0 / 6.0;
    %H_x(up) = - a * ( 11 * H(up) - 18 * H3(up) + 9 * H8(up) - 2 * H12(up) );
    %H_x(um) = a * ( 11 * H(um) - 18 * H2(um) + 9 * H7(um) - 2 * H11(um) );
    %H_y(vp) = - a * ( 11 * H(vp) - 18 * H4(vp) + 9 * H6(vp) - 2 * H10(vp) );
    %H_y(vm) = a * ( 11 * H(vm) - 18 * H1(vm) + 9 * H5(vm) - 2 * H9(vm) );
        

    % Maybe not neccessary for the second order?!
    %H_x(up & M==1) = - ( H(up & M==1) - H3(up & M==1) );
    %H_x(um & M==1) = H(um & M==1) - H2(um & M==1);
    %H_y(vp & M==1) = - ( H(vp & M==1) - H4(vp & M==1) );
    %H_y(vm & M==1) = H(vm & M==1) - H1(vm & M==1);

    % Boundary conditions applied on velocities.
    %M=zeros(ctr.imax,ctr.jmax);
    %M(1,:)=1;
    %M(ctr.imax,:)=1;
    %M(:,1)=1;
    %M(:,ctr.jmax)=1;
    
    %V1(M==1)=0;
    %V2(M==1)=0;
    %V3(M==1)=0;
    %V4(M==1)=0;

    % BOundary conditions.
    % CHECK BOUNDARY CONDITIONS.
    %u_stag(M==1) = 0.0;
    %v_stag(M==1) = 0.0;

    % Update ice thickness from velocity and surface gradients.
    %H_new = H + dtdx * ( u_stag .* H_x + v_stag .* H_y ) + R0;

    
    % Boundary conditions CalvingMIP.
    % CHECK BOUNDARY CONDITIONS
    %H_new(1,:) = H_new(2,:);
    %H_new(:,1) = H_new(:,2);
    %H_new(ctr.imax,:) = H_new(ctr.imax-1,:);
    %H_new(:,ctr.jmax) = H_new(:,ctr.jmax-1);

    %H_new(1,:) = H(1,:) + R0(1,:);
    %H_new(:,1) = H(:,1) + R0(:,1);
    %H_new(ctr.imax,:) = H(ctr.imax,:) + R0(ctr.imax,:);
    %H_new(:,ctr.jmax) = H(:,ctr.jmax) + R0(:,ctr.jmax);

    %H_new(1,:) = H(1,:);
    %H_new(:,1) = H(:,1);
    %H_new(ctr.imax,:) = H(ctr.imax,:);
    %H_new(:,ctr.jmax) = H(:,ctr.jmax);


    %M=zeros(ctr.imax,ctr.jmax);
    %M(1,:)=1;
    %M(ctr.imax,:)=1;
    %M(:,1)=1;
    %M(:,ctr.jmax)=1;

    qx = u .* H;
    qy = v .* H;

    a = 2:ctr.imax;
    b = 2:ctr.jmax;

    qy_y = qy - circshift(qy,[1 0]); % (i-1,j)
    qx_x = qx - circshift(qx,[0 1]); % (i,j-1)

    H_new(a,1) = H(a,2) - dtdx * ...
                            ( qx(a,3) - qx(a,2) + ...
                              qy_y(a,2) ) + R0(a,1);

    H_new(a,ctr.jmax) = H(a,ctr.jmax-1) - dtdx * ...
                            ( qx(a,ctr.jmax-1) - qx(a,ctr.jmax-2) + ...
                              qy_y(a,ctr.jmax-1) ) + R0(a,ctr.jmax);


    H_new(1,b) = H(2,b) - dtdx * ...
                            ( qx_x(2,b) + ...
                              qy(3,b) - qy(2,b) ) + R0(1,b);
  
    H_new(ctr.imax,b) = H(ctr.imax-1,b) - dtdx * ...
                            ( qx_x(ctr.imax-1,b) + ...
                              qy(ctr.imax-1,b) - qy(ctr.imax-2,b) ) + R0(ctr.imax,b);


%     for i=2:ctr.imax
%         
%         H_new(i,1) = H(i,2) - dtdx * ...
%                             ( qx(i,3) - qx(i,2) + ...
%                               qy(i,2) - qy(i-1,2) ) + R0(i,1);
% 
% 
%         H_new(i,ctr.jmax) = H(i,ctr.jmax-1) - dtdx * ...
%                                 ( qx(i,ctr.jmax-1) - qx(i,ctr.jmax-2) + ...
%                                   qy(i,ctr.jmax-1) - qy(i-1,ctr.jmax-1) ) + R0(i,ctr.jmax);
%     end
% 
% 
%     for j=2:ctr.jmax
%         
%         H_new(1,j) = H(2,j) - dtdx * ...
%                             ( qx(2,j) - qx(2,j-1) + ...
%                               qy(3,j) - qy(2,j) ) + R0(1,j);
% 
% 
%         H_new(ctr.imax,j) = H(ctr.imax-1,j) - dtdx * ...
%                                 ( qx(ctr.imax-1,j) - qx(ctr.imax-1,j-1) + ...
%                                   qy(ctr.imax-1,j) - qy(ctr.imax-2,j) ) + R0(ctr.imax,j); 
%     end

    % Corner is skipped in the loop.
    H_new(1,1) = H_new(1,2);



    % CHECK BOUNDARY CONDITIONS
    % EISMINT.
    %H_new(1,:) = 0.0;
    %H_new(:,1) = 0.0;
    %H_new(ctr.imax,:) = 0.0;
    %H_new(:,ctr.jmax) = 0.0;
 
    % Some relaxation avoids spurious results.
    rel = 0.25; % 0.25
    H = H_new * ( 1.0 - rel ) + rel * H;

    % Avoid spurious results.
    %H(H>3.0e3)=3.0e3;
    
    % Boundary conditions (identical to Frank's).
    if ctr.mismip>=1
        H(:,1) = H(:,3); % Divide.
        H(1,:) = H(3,:);

        if ctr.mismip==1
            H(ctr.imax,:) = H(1,:); % n-2,j - PBC at i=imax
        end
    end

end


