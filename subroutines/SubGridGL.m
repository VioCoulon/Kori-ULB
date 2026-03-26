function [betax, betay] = SubGridGL(beta2, H, ...
                                    HAF, MASK, glMASK, Hmx, Hmy, B, SLR, ctr, par)


    betax=0.5*(beta2+circshift(beta2,[0 -1]));
    betay=0.5*(beta2+circshift(beta2,[-1 0])); 


    HAF1 = circshift(HAF,[0 -1]); % (i,j+1)
    HAF2 = circshift(HAF,[-1 0]); % (i+1,j)

    % Working version.
    xm = (HAF>=0.0) & (HAF1<0.0);
    xp = (HAF<0.0) & (HAF1>=0.0);

    ym = (HAF>=0.0) & (HAF2<0.0);
    yp = (HAF<0.0) & (HAF2>=0.0);

    f_grnd_x = zeros(ctr.imax, ctr.jmax);
    f_grnd_y = zeros(ctr.imax, ctr.jmax);

    % Grounded fraction on u-grid.
    f_grnd_x(xm) = HAF(xm) ./ ( HAF(xm) - HAF1(xm) );
    f_grnd_y(ym) = HAF(ym) ./ ( HAF(ym) - HAF2(ym) );

    f_grnd_x(xp) = HAF1(xp) ./ ( HAF1(xp) - HAF(xp) );
    f_grnd_y(yp) = HAF2(yp) ./ ( HAF2(yp) - HAF(yp) );




    %f_grnd_x(xm) = HAF1(xm) ./ HAF(xm);
    %f_grnd_y(ym) = HAF2(ym) ./ HAF(ym);

    %f_grnd_x(xp) = HAF(xp) ./ HAF1(xp);
    %f_grnd_y(yp) = HAF(yp) ./ HAF2(yp);

    
    % Weighting from cells grounded fraction.
    %MASK1 = circshift(MASK, [0 -1]);
    %MASK2 = circshift(MASK, [0 1]);
    %MASK3 = circshift(MASK, [-1 0]);
    %MASK4 = circshift(MASK, [1 0]);

    %MASK6 = circshift(MASK, [1 -1]);
    %MASK7 = circshift(MASK, [-1 1]);
    %MASK8 = circshift(MASK, [-1 -1]);


    %M0 = (MASK == 1);
    %M1 = (MASK1 == 1);
    %M2 = (MASK2 == 1);
    %M3 = (MASK3 == 1);
    %M4 = (MASK4 == 1);
    %M5 = (MASK5 == 1);
    %M6 = (MASK6 == 1);
    %M7 = (MASK7 == 1);
    %M8 = (MASK8 == 1);


    %wx = ( MASK4 + MASK + MASK3 + MASK6 + MASK1 + MASK8 ) / 6.0;
    %ground = [MASK(i-1,j)==1, MASK(i,j)==1, MASK(i+1,j)==1, MASK(i-1,j+1)==1, MASK(i,j+1)==1, MASK(i+1,j+1)==1];


    %wy = ( MASK2 + MASK + MASK1 + MASK7 + MASK3 + MASK8 ) / 6.0;
    %ground = [MASK(i,j-1)==1, MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j-1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];

    %betax(xm) = betax(xm) ./ wx(xm);
    %betax(xp) = betax(xp) ./ wx(xp);

    %betay(xm) = betay(xm) ./ wy(xm);
    %betay(xp) = betay(xp) ./ wy(xp);



%     for i=1:ctr.imax
%        for j=1:ctr.jmax
%            %x-axis.
%            if (HAF(i,j)>=0.0) & (HAF1(i,j)<0.0)
% 
%                ground = [MASK(i-1,j)==1, MASK(i,j)==1, MASK(i+1,j)==1, MASK(i-1,j+1)==1, MASK(i,j+1)==1, MASK(i+1,j+1)==1];
%                %ground = [MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
% 
%                w = sum(ground) / 6.0;
%                
%                 
%                %betax(i,j) = w * beta2(i,j) + (1.0 - w) * beta2(i,j+1);
%                %betax(i,j) = beta2(i,j);
%                
%                ux(i,j+1) = 5*ux(i,j+1) / w;
%                ux(i,j) = 5*ux(i,j) / w;
% 
%                %ux(i,j-1) = 5*ux(i,j-1) / w;
% 
%            end
% 
% 
%            if (HAF(i,j)<0.0) & (HAF1(i,j)>=0.0)
% 
%                ground = [MASK(i-1,j)==1, MASK(i,j)==1, MASK(i+1,j)==1, MASK(i-1,j+1)==1, MASK(i,j+1)==1, MASK(i+1,j+1)==1];
%                %ground = [MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
% 
%                w = sum(ground) / 6.0;
% 
%                %betax(i,j) = w * beta2(i,j+1) + (1.0 - w) * beta2(i,j);
%                %betax(i,j) = beta2(i,j+1);
% 
%                ux(i,j-1) = 5*ux(i,j-1) / w;
%                ux(i,j) = 5*ux(i,j) / w;
% 
%                %ux(i,j+1) = 5*ux(i,j+1) / w;
% 
%            end
% 
% 
% 
%            %y-axis.
%            if (HAF(i,j)>=0.0) & (HAF2(i,j)<0.0)
% 
%                ground = [MASK(i,j-1)==1, MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j-1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
%                %ground = [MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
% 
%                w = sum(ground) / 6.0;
% 
%                %betay(i,j) = w * beta2(i,j) + (1.0 - w) * beta2(i+1,j);
%                %betay(i,j) = beta2(i,j);
%                uy(i+1,j) = 5* uy(i+1,j) / w;
%                uy(i,j) = 5*uy(i,j) / w;
% 
%                %uy(i-1,j) = 5*uy(i-1,j) / w;
% 
%            end
% 
% 
%            if (HAF(i,j)<0.0) & (HAF2(i,j)>=0.0)
% 
%                ground = [MASK(i,j-1)==1, MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j-1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
%                %ground = [MASK(i,j)==1, MASK(i,j+1)==1, MASK(i+1,j)==1, MASK(i+1,j+1)==1];
% 
%                w = sum(ground) / 6.0;
%                
% 
%                %betay(i,j) = w * beta2(i+1,j) + (1.0 - w) * beta2(i,j);
%                %betay(i,j) = beta2(i+1,j);
% 
%                uy(i-1,j) = 5*uy(i-1,j) / w;
%                uy(i,j) = 5*uy(i,j) / w;
% 
%                %uy(i+1,j) = 5*uy(i+1,j) / w;
% 
%            end
% 
%        end 
%     end


    



    % Stagger beta for interpolation. 
    %beta2_x = circshift(beta2,[0 -1]);  % (i,j+1)
    %beta2_y = circshift(beta2,[-1 0]);  % (i+1,j)

    %betax_x = circshift(betax,[0 -1]);  % (i,j+1)
    %betay_y = circshift(betay,[-1 0]);


    % Weighting term as a function of grounded fraction 
    wt_x = f_grnd_x;
    wt_y = f_grnd_y;

    % Full advance.
    % Friction evaluated in u-grid from interpolation.
    %betax(xm) = wt_x(xm) .* beta2_x(xm) + (1.0-wt_x(xm)) .* beta2(xm);
    %betax(xp) = (1.0-wt_x(xp)) .* beta2_x(xp) + wt_x(xp) .* beta2(xp);

    %betay(ym) = wt_y(ym) .* beta2_y(ym) + (1.0-wt_y(ym)) .* beta2(ym);
    %betay(yp) = (1.0-wt_y(yp)) .* beta2_y(yp) + wt_y(yp) .* beta2(yp);

    %xm2 = circshift(xm, [0 -1]); % (HAF>=0.0) & (HAF1<0.0);
    %xp2 = circshift(xp, [0 1]);% (HAF<0.0) & (HAF1>=0.0);

    %ym2 = circshift(ym, [-1 0]); %(HAF>=0.0) & (HAF2<0.0);
    %yp2 = circshift(yp, [1 0]); %(HAF<0.0) & (HAF2>=0.0);


    % GL advance is limited. SOMETHING WRONG???
    %betax(xm) = wt_x(xm) .* beta2(xm) + (1.0-wt_x(xm)) .* beta2_x(xm);
    %betax(xp) = (1.0-wt_x(xp)) .* beta2(xp) + wt_x(xp) .* beta2_x(xp);

    %betay(ym) = wt_y(ym) .* beta2(ym) + (1.0-wt_y(ym)) .* beta2_y(ym);
    %betay(yp) = (1.0-wt_y(yp)) .* beta2(yp) + wt_y(yp) .* beta2_y(yp);

    
    % BEST ONE SO FAR.
    betax(xm) = wt_x(xm) .* betax(xm);
    betax(xp) = (1.0-wt_x(xp)) .* betax(xp);

    betay(ym) = wt_y(ym) .* betay(ym);
    betay(yp) = (1.0-wt_y(yp)) .* betay(yp);


    %betax(xm) = wt_x(xm) .* beta2(xm);
    %betax(xp) = wt_x(xp) .* beta2_x(xp);

    %betay(ym) = wt_y(ym) .* beta2(ym);
    %betay(yp) = wt_y(yp) .* beta2_y(yp);


    % Try correcting the velicties that will be used to calculate ice thickness.
    % Similar to what schoof does.
    %ux_x = circshift(ux,[0 -1]);  % (i,j+1)
    %uy_y = circshift(uy,[-1 0]);

    %ux(xm) = (1.0+wt_x(xm)) .* ux(xm);
    %ux(xp) = (1.0-wt_x(xp)) .* ux(xp);

    %uy(ym) = (1.0+wt_y(ym)) .* uy(ym) ;
    %uy(yp) = (1.0+wt_y(yp)) .* uy(yp);

    %max(wt_x(xm),[], 'All')
    %min(wt_x(xm),[], 'All')


    %xm2 = circshift(xm, [0 1]); % (HAF>=0.0) & (HAF1<0.0);
    %xp2 = circshift(xp, [0 -1]);% (HAF<0.0) & (HAF1>=0.0);

    %ym2 = circshift(ym, [1 0]); %(HAF>=0.0) & (HAF2<0.0);
    %yp2 = circshift(yp, [-1 0]); %(HAF<0.0) & (HAF2>=0.0);


    %ux(xm2) = ux(xm2) ./ wt_x(xm);
    %ux(xp2) = ux(xp2) ./ wt_x(xp);

    %uy(ym2) = uy(ym2) ./ wt_y(ym);
    %uy(yp2) = uy(yp2) ./ wt_y(yp);


    %betax = 1e3 * betax;
    %betay = 1e3 * betay;


    % Use structure from schoofing. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conditions for GL in x-direction and ux>0
%     Mgl=zeros(ctr.imax,ctr.jmax);
% 
%     glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
%     glMASK0=circshift(glMASK,[0 1]); % glMASK(i,j-1), upstream
%     ux1=circshift(ux,[0 1]); % ux(i,j-1), upstream velocity
%     HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
%     B1=circshift(B,[0 -1]); % B(i,j+1)
% 
%     Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & ux>0 & ux1>0)=1; % u-grid(i,j) 
%     fracgx=HAF./(HAF-HAF1);
%     hbgx=(1.-fracgx).*B+fracgx.*B1;
%     hgx=(SLR-hbgx)*par.rhow/par.rho;
%     Mgl(hgx<0)=0;
% 
% 
%     % ice shelf grid cell downstream of GL
%     Mgl=circshift(Mgl,[0 1]); % condition for i,j+1
%     Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmx==0)=0; % only apply if shelf exists
% 
%     wt = abs(fracgx);
% 
%     Mgl=circshift(Mgl,[0 1]); % condition for i,j+1
%     Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmx==0)=0; % only apply if shelf exists
%     ux(Mgl==1) = ux(Mgl==1) ./ wt(Mgl==1);    
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % conditions for GL in x-direction and ux<0
% 
%     Mgl=zeros(ctr.imax,ctr.jmax);
% 
%     glMASK1=circshift(glMASK,[0 -1]); % glMASK(i,j+1)
%     glMASK0=circshift(glMASK,[0 -2]); % glMASK(i,j+2), upstream
%     ux1=circshift(ux,[0 -1]);   % ux(i,j+1), upstream velocity
%     HAF1=circshift(HAF,[0 -1]); % HAF(i,j+1)
%     B1=circshift(B,[0 -1]); % B(i,j+1)
%     Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & ux<0 & ux1<0)=1; % u-grid(i,j)
%     fracgx=HAF1./(HAF1-HAF);
%     hbgx=(1.-fracgx).*B1+fracgx.*B;
%     hgx=(SLR-hbgx)*par.rhow/par.rho;
%     Mgl(hgx<0)=0;
%     %VL: GL in x-direction, u<0 --> jh=jGL, jg=jh+1, ih=ig=iGL
%     
%     % ice shelf grid cell downstream of GL
%     Mgl=circshift(Mgl,[0 -1]); % condition for i,j-1
%     Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmx==0)=0; % only apply if shelf exists
%     
% 
%     wt = abs(fracgx);
% 
%     Mgl=circshift(Mgl,[1 0]); % condition for i+1,j
%     Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmy==0)=0; % only apply if shelf exists
% 
% 
%     ux(Mgl==1) = ux(Mgl==1) ./ wt(Mgl==1);  
% 
% 
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % conditions for GL in y-direction and uy>0
%    
%     Mgl=zeros(ctr.imax,ctr.jmax);
%     
%     glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
%     glMASK0=circshift(glMASK,[1 0]); % glMASK(i-1,j), upstream
%     uy1=circshift(uy,[1 0]);    % uy(i-1,j), upstream velocity
%     HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
%     B1=circshift(B,[-1 0]); % B(i+1,j)
%     
%     Mgl(glMASK==2 & glMASK1>=3 & glMASK0<=2 & uy>0 & uy1>0)=1; % v-grid(i,j)
%     fracgy=HAF./(HAF-HAF1);
%     hbgy=(1.-fracgy).*B+fracgy.*B1;
%     hgy=(SLR-hbgy)*par.rhow/par.rho;
%     Mgl(hgy<0)=0;
% 
%     %uy(Mgl==1) = uy(Mgl==1) ./ wt(Mgl==1);
%     
%     
%     % ice shelf grid cell downstream of GL
%     Mgl=circshift(Mgl,[1 0]); % condition for i+1,j
%     Mgl(glMASK1<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmy==0)=0; % only apply if shelf exists
% 
%     wt = abs(fracgy);
% 
%     uy(Mgl==1) = uy(Mgl==1) ./ wt(Mgl==1); 
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
%     % conditions for GL in y-direction and uy<0
%     Mgl=zeros(ctr.imax,ctr.jmax);
%     
%     glMASK1=circshift(glMASK,[-1 0]); % glMASK(i+1,j)
%     glMASK0=circshift(glMASK,[-2 0]); % glMASK(i+2,j), upstream
%     uy1=circshift(uy,[-1 0]);   % uy(i+1,j), upstream velocity
%     HAF1=circshift(HAF,[-1 0]); % HAF(i+1,j)
%     B1=circshift(B,[-1 0]); % B(i+1,j)
%     Mgl(glMASK>=3 & glMASK1==2 & glMASK0<=2 & uy<0 & uy1<0)=1; % v-grid(i,j)
%     fracgy=HAF1./(HAF1-HAF);
%     hbgy=(1.-fracgy).*B1+fracgy.*B;
%     hgy=(SLR-hbgy)*par.rhow/par.rho;
%     Mgl(hgy<0)=0;
%     
%     % ice shelf grid cell downstream o3
%     Mgl=circshift(Mgl,[-1 0]); % condition for i-1,j
%     Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmy==0)=0; % only apply if shelf exists
% 
% 
%     wt = abs(fracgy);
% 
%     Mgl=circshift(Mgl,[-1 0]); % condition for i-1,j
%     Mgl(glMASK<=2)=0; % only apply if downstream grid cell is floating
%     Mgl(Hmy==0)=0; % only apply if shelf exists
% 
% 
%     uy(Mgl==1) = uy(Mgl==1) ./ wt(Mgl==1); 


end
