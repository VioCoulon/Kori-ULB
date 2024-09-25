function H=SolverIceThickness_optimised(Mb,H,u,v,glMASK,ctr)

% Kori-ULB
% Optimized explicit solver of the ice thickness equation

    % Daniel's explicit calculation of ice thickness.
    %dt_CFL = 0.75 * ctr.delta * ( (1.0./max(u)) + (1.0/max(v)) );
    %dt = min(ctr.dt,dt_CFL);
    dtdx=ctr.dt/ctr.delta;

    R0 = Mb*ctr.dt;
    % This cannot be according to CalvMIP Exp 1!
    %R0(MASK==0)=0.0;
    %R0=zeros(ctr.imax,ctr.jmax);

    % Staggered definitions.
    H1=circshift(H,[-1 0]); % (i+1,j)
    H2=circshift(H,[0 -1]); % (i,j+1)
    H3=circshift(H,[0 1]); % (i,j-1)
    H4=circshift(H,[1 0]); % (i-1,j)

    %H5=circshift(H,[-2 0]); % (i+2,j)
    %H6=circshift(H,[2 0]); % (i-2,j)
    %H7=circshift(H,[0 -2]); % (i,j+2)
    %H8=circshift(H,[0 2]); % (i,j-2)
    
    %H9=circshift(H,[-3 0]); % (i+3,j)
    %H10=circshift(H,[3 0]); % (i-3,j)
    %H11=circshift(H,[0 -3]); % (i,j+3)
    %H12=circshift(H,[0 3]); % (i,j-3)

    
    %H_x = zeros(ctr.imax,ctr.jmax);
    %H_y = zeros(ctr.imax,ctr.jmax);

    % Stagger velocities.
    %u1 = circshift(u,[0 1]); % (i,j-1)
    %v1 = circshift(v,[1 0]); % (i-1,j)
    %u_stag = 0.5 * ( u + u1 );
    %v_stag = 0.5 * ( v + v1 );

    % No staggering if using the same grid.
    u_stag = u;
    v_stag = v;

    
    % Boolean indexes.
    up = u_stag>0.0;
    um = u_stag<0.0;
    vp = v_stag>0.0;
    vm = v_stag<0.0;
    calv = glMASK==5;


    % ORIGINAL FORMULATION. First prder derivatives.
    %H_x(up) = - ( H(up) - H3(up) );
    %H_x(um) = H(um) - H2(um);
    %H_y(vp) = - ( H(vp) - H4(vp) );
    %H_y(vm) = H(vm) - H1(vm);

    % IT WORKS NOW!
    %H_x(up) = - ( H2(up) - H(up) );
    %H_x(um) = H(um) - H2(um);
    %H_y(vp) = - ( H1(vp) - H(vp) );
    %H_y(vm) = H(vm) - H1(vm);

    % Working.
    H_x = H - H2;
    H_y = H - H1;
    
    %H_x = H3 - H;
    %H_y = H4 - H;

    %H_x(up) = - 0.5 * ( H2(up) - H3(up) );
    %H_x(um) = 0.5 * ( H3(um) - H2(um) );
    %H_y(vp) = - 0.5 * ( H1(vp) - H4(vp) );
    %H_y(vm) = 0.5 * ( H4(vm) - H1(vm) );

    %H_x(up) = - 0.5 * ( H(up) - 4.0 * H3(up) + 3.0 * H8(up));
    %H_x(um) = 0.5 * ( H(um) - 4.0 * H2(um) + 3.0 * H7(um) );
    %H_y(vp) = - 0.5 * ( H(vp) - 4.0 * H4(vp) + 3.0 * H6(vp));
    %H_y(vm) = 0.5 * ( H(vm) - 4.0 * H1(vm) + 3.0 * H5(vm) );



    %H_x(up) = - ( H2(up) - H(up) );
    %H_x(um) =  H3(um) - H(um);
    %H_y(vp) = - ( H1(vp) - H(vp) );
    %H_y(vm) = H4(vm) - H(vm);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use fluxes instead.
    %q_x = zeros(ctr.imax,ctr.jmax);
    %q_y = zeros(ctr.imax,ctr.jmax);

    %q_x(up) = u_stag(up) .* H_x(up);
    %q_x(um) = u_stag(um) .* H_x(um); 
    %q_y(vp) = v_stag(vp) .* H_y(vp);
    %q_y(vm) = v_stag(vm) .* H_y(vm);

    q_x = u_stag .* H_x;
    q_y = v_stag .* H_y;

    % Update ice thickness from velocity and surface gradients.
    qx_stag = 0.5 * ( q_x + circshift(q_x, [0 1]) );
    qy_stag = 0.5 * ( q_y + circshift(q_y, [1 0]) );

    
    H_new = H + dtdx * ( qx_stag + qy_stag ) + R0;
    %H_new = H + dtdx * ( q_x + q_y ) + R0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calving front.
    %H1=circshift(H,[-1 0]); % (i+1,j)
    %H2=circshift(H,[0 -1]); % (i,j+1)
    %H3=circshift(H,[0 1]); % (i,j-1)
    %H4=circshift(H,[1 0]); % (i-1,j)

    u1 = circshift(u_stag,[0 1]); % (i,j-1)
    %H_new(calv & up) = H(calv & up) + dtdx * ...
    %                        ( u1(calv & up) .* ( H3(calv & up) - H(calv & up) ) +...
    %                            q_y(calv & up) ) +...
    %                        R0(calv & up);

    % CHANGE THE FLUX IN Y DIRECTION HERE!!!!!
    H_new(calv & up) = H(calv & up) + dtdx * ...
                                ( u1(calv & up) .* ( H3(calv & up) - H(calv & up) ) +...
                                    qy_stag(calv & up) ) +...
                                R0(calv & up);


    H_new(calv & um) = H(calv & um) + dtdx * ...
                            ( u_stag(calv & um) .* ( H(calv & um) - H2(calv & um) ) +...
                                qy_stag(calv & um) ) +...
                            R0(calv & um);

    v1 = circshift(v_stag,[1 0]); % (i-1,j)
    H_new(calv & vp) = H(calv & vp) + dtdx * ...
                            ( qx_stag(calv & vp) + ...
                                v1(calv & vp) .* ( H4(calv & vp) - H(calv & vp) ) ) +...
                            R0(calv & vp);

    H_new(calv & vm) = H(calv & vm) + dtdx * ...
                            ( qx_stag(calv & vm) + ...
                                v_stag(calv & vm) .* ( H(calv & vm) - H1(calv & vm) ) ) +...
                            R0(calv & vm);
                                 
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

    
    % Boundary conditions CalvingMIP.ç
    % CHECK BOUNDARY CONDITIONS
    H_new(1,:) = H_new(2,:);
    H_new(:,1) = H_new(:,2);
    H_new(ctr.imax,:) = H_new(ctr.imax-1,:);
    H_new(:,ctr.jmax) = H_new(:,ctr.jmax-1);

    %H_new(1,:) = H(1,:) + R0(1,:);
    %H_new(:,1) = H(:,1) + R0(:,1);
    %H_new(ctr.imax,:) = H(ctr.imax,:) + R0(ctr.imax,:);
    %H_new(:,ctr.jmax) = H(:,ctr.jmax) + R0(:,ctr.jmax);

    % CHECK BOUNDARY CONDITIONS
    % EISMINT.
    %H(1,:) = 0.0;
    %H(:,1) = 0.0;
    %H(ctr.imax,:) = 0.0;
    %H(:,ctr.jmax) = 0.0;
 
    % Some relaxation avoids spurious results.
    rel = 0.25; % 0.25
    H = H_new * ( 1.0 - rel ) + rel * H;


    % Avoid spurious results.
    H(H>3.0e3)=3.0e3;
    
    % Boundary conditions (identical to Frank's).
    if ctr.mismip>=1
        H(:,1) = H(:,3); % Divide.
        H(1,:) = H(3,:);

        if ctr.mismip==1
            H(ctr.imax,:) = H(1,:); % n-2,j - PBC at i=imax
        end
    end

end


