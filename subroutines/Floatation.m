function [HAF,MASK,HB,sn]=Floatation(par,B,SLR,H,MASK)

% Kori-ULB
% Determination of MASK and height of the bottom of ice shelves based on
% floatation. Height above buoyancy is also returned.

    % Artificially induced GL retreat.
    epsilon=10.0;

    H(MASK==0 & H<=par.SeaIceThickness)=0;
    HAF=B-SLR+H*par.rho/par.rhow;
    MASK(HAF<0)=0;   % MASK(HAF<0)=0;
    MASK(HAF>=0)=1;  % MASK(HAF>=0)=1;
    %MASK(HAF<=epsilon)=0; 
    %MASK(HAF>epsilon)=1;  % MASK(HAF>=0)=1;
    HB=max(SLR-par.rho/par.rhow*H,B);
    sn=HB+H;


    % Daniel.
    % Correction to allow for retreat.
    %B1   = circshift(B,[0 1]);
    %H2   = circshift(H,[0 -1]);
    %HAF1 = circshift(HAF,[0 1]);
    %HAF2 = circshift(HAF,[0 -1]);
    %HAF12=B1-SLR+H2*par.rho/par.rhow;

    %a12 = (HAF1>=0) & (HAF2<0);

    %MASK((HAF12<0) & (a12))=0;   % MASK(HAF<0)=0;
    %MASK((HAF12>=0) & (a12))=1;  % MASK(HAF>=0)=1;



    %B3   = circshift(B,[1 0]);
    %H4   = circshift(H,[-1 0]);
    %HAF3 = circshift(HAF,[1 0]);
    %HAF4 = circshift(HAF,[-1 0]);
    %HAF34=B3-SLR+H4*par.rho/par.rhow;

    %a34 = (HAF3>=0) & (HAF4<0);

    %MASK((HAF34<0) & (a34))=0;   % MASK(HAF<0)=0;
    %MASK((HAF34>=0) & (a34))=1;  % MASK(HAF>=0)=1;



    %B5   = circshift(B,[0 -1]);
    %H6   = circshift(H,[0 1]);
    %HAF5 = circshift(HAF,[0 -1]);
    %HAF6 = circshift(HAF,[0 1]);
    %HAF56=B5-SLR+H6*par.rho/par.rhow;
    
    %a56 = (HAF5>=0) & (HAF6<0);

    %MASK((HAF56<0) & (a56))=0;   % MASK(HAF<0)=0;
    %MASK((HAF56>=0) & (a56))=1;  % MASK(HAF>=0)=1;



    %B7   = circshift(B,[-1 0]);
    %H8   = circshift(H,[1 0]);
    %HAF7 = circshift(HAF,[-1 0]);
    %HAF8 = circshift(HAF,[1 0]);
    %HAF78=B7-SLR+H8*par.rho/par.rhow;

    %a78 = (HAF7>=0) & (HAF8<0);

    %MASK((HAF78<0) & (a78))=0;   % MASK(HAF<0)=0;
    %MASK((HAF78>=0) & (a78))=1;  % MASK(HAF>=0)=1;
    
end