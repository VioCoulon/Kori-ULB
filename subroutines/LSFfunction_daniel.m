function LSF=LSFfunction_daniel(LSF,ctr,u,v,node,nodes,VM,MASK,glMASK,X,Y,LSFo)

% Kori-ULB
% Calculate the Level Set Function (LSF) for following the calving front.
% Used in the calving algorihms
% Still under development


    % Daniel's explicit calculation of LSF.
    dtdx=ctr.dt/ctr.delta;

    R0 = LSF;
    %LSF_now = LSF;

    LSF1=circshift(LSF,[-1 0]); % (i+1,j)
    LSF2=circshift(LSF,[0 -1]); % (i,j+1)
    %LSF3=circshift(LSF,[0 1]); % (i,j-1)
    %LSF4=circshift(LSF,[1 0]); % (i-1,j)

    %LSF_x = 0.5 * (LSF3 - LSF2);
    %LSF_y = 0.5 * (LSF4 - LSF1);

    LSF_x = 0.5 * ( LSF + LSF2 );
    LSF_y = 0.5 * ( LSF + LSF1 );
    
    q_x = u .* LSF_x;
    q_y = v .* LSF_y;
    
    q_x1 = circshift(q_x,[0 1]); % (i,j-1)
    q_y1 = circshift(q_y,[1 0]); % (i,j-1)

    q_x = q_x - q_x1;
    q_y = q_y - q_y1;

    % WRONG??? U AND V ALSO CHANGE SPATIALLY.
    % CHANGE THIS TO ACCOUNT FOR SPATIALLY VARIABLE VELOCITY!!!

    %LSF = LSF + dtdx * ( u .* LSF_x + v .* LSF_y );
    LSF = LSF - dtdx * ( q_x + q_y );

    % Apply bounds.
    LSF(LSF<-1.0) = -1.0;
    LSF(LSF>1.0)  = 1.0;

    % TEST ON CALVING FRONT APPROACHING GL.
    %M = zeros(ctr.imax,ctr.jmax);

    M1 = circshift(MASK,[3 3]);
    M2 = circshift(MASK,[3 -3]);
    M3 = circshift(MASK,[-3 3]);
    M4 = circshift(MASK,[-3 -3]);

    %M5 = circshift(MASK,[0 2]);
    %M6 = circshift(MASK,[2 0]);
    %M7 = circshift(MASK,[0 -2]);
    %M8 = circshift(MASK,[-2 0]);



    %a = (MASK==1)|(M1==1)|(M2==1)|(M3==1)|(M4==1)...
    %             |(M5==1)|(M6==1)|(M7==1)|(M8==1);

    a = (MASK==1)|(M1==1)|(M2==1)|(M3==1)|(M4==1);
    M(a) = 1;


    % Daniel: calving front cannot retreat further than the GL by definition.
    %LSF(MASK==1) = R0(MASK==1);
    LSF(M==1) = R0(M==1);

    % Calving front cannot advance more than initial position.
    LSF(LSFo<0)=LSFo(LSFo<0);

end

