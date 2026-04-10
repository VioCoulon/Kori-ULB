clear all
close all
cd 'C:\Users\User\Desktop\ULB\CaMIP-Exp1'

ExpName='CalvingMIP-Exp1-Kori.nc';
 
if exist(ExpName,'file')==2
  delete(ExpName);
end

 % % 
%   for t = 0:1:999
%     N='Thule-Alg-2-C_';
%     T=sprintf('%03d', t);
%     NT=strcat(N,T);
%     load(NT)
load Exp1_toto.mat

Time=0;
 X=-795e3:10e3:795e3;
 Y=-795e3:10e3:795e3;

t=0;
VX(:,:,t+1)=ux;
VY(:,:,t+1)=uy;
H(:,:,t+1)=H;
lsf(:,:,t+1)=LSF;

MASK(MASK==0)=MASK(MASK==0)+1;
MASK(LSF<0)=3;
Mask(:,:,t+1)=MASK;

[X1,Y1,P1]=improfile(H,[80.5 80.5],[160 2],'bilinear'  );
[X2,Y2,P2]=improfile(H,[1 160],[160 1],'bilinear' );
[X3,Y3,P3]=improfile(H,[160 2] ,[80.5 80.5],'bilinear' );
[X4,Y4,P4]=improfile(H,[1 160],[1 160],'bilinear' );

[R1]=improfile(ux,[80.5 80.5],[160 2],'bilinear'  );
[R2]=improfile(ux,[1 160],[160 1],'bilinear' );
[R3]=improfile(ux,[160 2] ,[80.5 80.5],'bilinear' );
[R4]=improfile(ux,[1 160],[1 160],'bilinear' );

[S1]=improfile(uy,[80.5 80.5],[160 2],'bilinear'  );
[S2]=improfile(uy,[1 160],[160 1],'bilinear' );
[S3]=improfile(uy,[160 2] ,[80.5 80.5],'bilinear' );
[S4]=improfile(uy,[1 160],[1 160],'bilinear' );

nccreate(ExpName,'VX','Dimensions',{'X' 160 'Y' 160 'VX (ma-1)' numel(Time)})
nccreate(ExpName,'VY','Dimensions',{'X' 160 'Y' 160 'VY (ma-1)' numel(Time)})
nccreate(ExpName,'H','Dimensions',{'X' 160 'Y' 160 'H (m)' numel(Time)})
nccreate(ExpName,'Mask','Dimensions',{'X' 160 'Y' 160 'Ice mask' numel(Time)})
% nccreate(ExpName,'LSF','Dimensions',{'X' 160 'Y' 160 'Time (a)' numel(Time)})

nccreate(ExpName,'Time','Dimensions',{'Time (a)' numel(Time)})
nccreate(ExpName,'X','Dimensions',{'X position (m)' 160})
nccreate(ExpName,'Y','Dimensions',{'Y position (m)' 160})

nccreate(ExpName,'Profile_A_Thickness','Dimensions',{'Profile A ice thickness (m)' numel(P1)})
nccreate(ExpName,'Profile_A_X','Dimensions',{'Profile A point X coord (m)' numel(P1)})
nccreate(ExpName,'Profile_A_Y','Dimensions',{'Profile A point Y coord' numel(P1)})
nccreate(ExpName,'Profile_A_UVel','Dimensions',{'Profile A UVel (m a-1)' numel(P1)})
nccreate(ExpName,'Profile_A_VVel','Dimensions',{'Profile A VVel (m a-1)' numel(P1)})

nccreate(ExpName,'Profile_B_Thickness','Dimensions',{'Profile B ice thickness (m)' numel(P2)})
nccreate(ExpName,'Profile_B_X','Dimensions',{'Profile B point X coord (m)' numel(P2)})
nccreate(ExpName,'Profile_B_Y','Dimensions',{'Profile B point Y coord (m)' numel(P2)})
nccreate(ExpName,'Profile_B_UVel','Dimensions',{'Profile B UVel (m a-1)' numel(P2)})
nccreate(ExpName,'Profile_B_VVel','Dimensions',{'Profile B VVel (m a-1)' numel(P2)})

nccreate(ExpName,'Profile_C_Thickness','Dimensions',{'Profile C ice thickness (m)' numel(P3)})
nccreate(ExpName,'Profile_C_X','Dimensions',{'Profile C point X coord (m)' numel(P3)})
nccreate(ExpName,'Profile_C_Y','Dimensions',{'Profile C point Y coord (m)' numel(P3)})
nccreate(ExpName,'Profile_C_UVel','Dimensions',{'Profile C UVel (m a-1)' numel(P3)})
nccreate(ExpName,'Profile_C_VVel','Dimensions',{'Profile C VVel (m a-1)' numel(P3)})

nccreate(ExpName,'Profile_D_Thickness','Dimensions',{'Profile D ice thickness (m)' numel(P4)})
nccreate(ExpName,'Profile_D_X','Dimensions',{'Profile D point X coord (m)' numel(P4)})
nccreate(ExpName,'Profile_D_Y','Dimensions',{'Profile D point Y coord (m)' numel(P4)})
nccreate(ExpName,'Profile_D_UVel','Dimensions',{'Profile D UVel (m a-1)' numel(P4)})
nccreate(ExpName,'Profile_D_VVel','Dimensions',{'Profile D VVel (m a-1)' numel(P4)})

ncwrite(ExpName,'Time',Time)
ncwrite(ExpName,'X',X)
ncwrite(ExpName,'Y',Y)
ncwrite(ExpName,'VX',VX)
ncwrite(ExpName,'VY',VY)
ncwrite(ExpName,'H',H)
ncwrite(ExpName,'Mask',Mask)
% ncwrite(ExpName,'LSF',lsf)

ncwrite(ExpName,'Profile_A_Thickness',P1)
ncwrite(ExpName,'Profile_A_X',X1)
ncwrite(ExpName,'Profile_A_Y',Y1)
ncwrite(ExpName,'Profile_A_Thickness',R1)
ncwrite(ExpName,'Profile_A_Thickness',S1)

ncwrite(ExpName,'Profile_B_Thickness',P2)
ncwrite(ExpName,'Profile_B_X',X2)
ncwrite(ExpName,'Profile_B_Y',Y2)
ncwrite(ExpName,'Profile_B_Thickness',R2)
ncwrite(ExpName,'Profile_B_Thickness',S2)

ncwrite(ExpName,'Profile_C_Thickness',P3)
ncwrite(ExpName,'Profile_C_X',X3)
ncwrite(ExpName,'Profile_C_Y',Y3)
ncwrite(ExpName,'Profile_C_Thickness',R3)
ncwrite(ExpName,'Profile_C_Thickness',S3)

ncwrite(ExpName,'Profile_D_Thickness',P4)
ncwrite(ExpName,'Profile_D_X',X4)
ncwrite(ExpName,'Profile_D_Y',X4)
ncwrite(ExpName,'Profile_D_Thickness',R4)
ncwrite(ExpName,'Profile_D_Thickness',S4)

ncdisp(ExpName)






