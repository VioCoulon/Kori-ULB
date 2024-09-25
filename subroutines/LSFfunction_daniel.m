function LSF=LSFfunction(LSF,ctr,u,v,node,nodes,VM,MASK,glMASK,X,Y,LSFo)

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
    LSF3=circshift(LSF,[0 1]); % (i,j-1)
    LSF4=circshift(LSF,[1 0]); % (i-1,j)

    LSF_x = 0.5 * (LSF3 - LSF2);
    LSF_y = 0.5 * (LSF4 - LSF1);

    LSF = LSF + dtdx * ( u .* LSF_x + v .* LSF_y );

    % Apply bounds.
    LSF(LSF<-1.0) = -1.0;
    LSF(LSF>1.0)  = 1.0;

    % Daniel: calving front cannot retreat further than the GL by definition.
    LSF(MASK==1) = R0(MASK==1);

end


