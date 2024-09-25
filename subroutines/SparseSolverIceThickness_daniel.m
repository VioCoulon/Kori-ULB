function H=SparseSolverIceThickness_daniel(Mb,H,dtdx,u,v,ctr)

% Kori-ULB
% Optimized explicit solver of the ice thickness equation

    % Daniel's explicit calculation of ice thickness.
    %dt_CFL = 0.75 * ctr.delta * ( (1.0./max(u)) + (1.0/max(v)) );
    %dt = min(ctr.dt,dt_CFL);
    %dtdx=dt/ctr.delta;



    R0 = Mb*ctr.dt;
    % This cannot be according to CalvMIP Exp 1!
    %R0(MASK==0)=0.0;
    %R0=zeros(ctr.imax,ctr.jmax);


            
%     H_new = H;
%     for i=2:ctr.imax-1
%         for j=2:ctr.jmax-1
%             if u_stag(i,j) > 0.0 % u_stag(i,j) > 0.0
%                 
%                 % It works.
%                 %H_x = -(H(i,j+1) - H(i,j));
% 
%                 %Symmetric attempts..
%                 %H_x = -0.5*(H(i,j+1) - H(i,j-1));
% 
%                 % It works.
%                 H_x = -(H(i,j) - H(i,j-1));
% 
%             end
%             if u_stag(i,j) < 0.0
% 
%                 % It works.
%                 %H_x = H(i,j) - H(i,j+1);
% 
%                 %Symmetric attempts..
%                 %H_x = 0.5*(H(i,j-1) - H(i,j+1));
% 
%                 % It works.
%                 H_x = H(i,j) - H(i,j+1);
% 
%             end
%             if v_stag(i,j) > 0.0
% 
%                 %Works.
%                 %H_y = -(H(i+1,j) - H(i,j));
% 
%                 %Symmetric attempts..
%                 %H_y = -0.5*(H(i+1,j) - H(i-1,j));
% 
%                 %Works.
%                 H_y = -(H(i,j) - H(i-1,j));
% 
%             end
%             if v_stag(i,j) < 0.0
% 
%                 % It wkrs.
%                 %H_y = H(i,j) - H(i+1,j);
% 
%                 %Symmetric attempts..
%                 %H_y = 0.5*(H(i-1,j) - H(i+1,j));
% 
%                 % It wkrs.
%                 H_y = H(i,j) - H(i+1,j);
% 
%             end
% 
%             % Works.
%             H_new(i,j) = H(i,j) + dtdx * ( u_stag(i,j)*H_x + v_stag(i,j)*H_y ) + R0(i,j);
% 
% 
%         end 
%     end


    % Staggered definitions.
    H1=circshift(H,[-1 0]); % (i+1,j)
    H2=circshift(H,[0 -1]); % (i,j+1)
    H3=circshift(H,[0 1]); % (i,j-1)
    H4=circshift(H,[1 0]); % (i-1,j)

    u1 = circshift(u,[0 1]); % (i,j-1)
    v1 = circshift(v,[1 0]); % (i-1,j)

    H_x = zeros(ctr.imax,ctr.jmax);
    H_y = zeros(ctr.imax,ctr.jmax);

    % Fractor 0.5 * (u + u1)
    u_stag = 0.5 * ( u + u1 );
    v_stag = 0.5 * ( v + v1 );

    % NO Factor 0.5!
    up = u_stag>0.0;
    um = u_stag<0.0;
    vp = v_stag>0.0;
    vm = v_stag<0.0;

    H_x(up) = - ( H(up) - H3(up) );
    H_x(um) = H(um) - H2(um);
    H_y(vp) = - ( H(vp) - H4(vp) );
    H_y(vm) = H(vm) - H1(vm);

    % Factor 0.25 from 0.5 of H_x and 0.5 from u_stag.
    % 0.5 for dt.
    H_new = H + dtdx * ( u_stag.*H_x + v_stag.*H_y ) + R0;
 
    % Some relaxation.
    rel = 0.2; % 0.25
    H = H_new * ( 1.0 - rel ) + rel * H;

    % Boundary conditions CalvingMIP.
    H(1,:) = H(2,:);
    H(:,1) = H(:,2);
    H(ctr.imax,:) = H(ctr.imax-1,:);
    H(:,ctr.jmax) = H(:,ctr.jmax-1);


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


