function Thule


parent_path = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/';
%path_in     = [parent_path, 'ice_data/calvingMIP/Exp3/dx_2km/'];
path_in     = [parent_path, 'output/calvingMIP/Exp3-4/dx_2km/OceanVisc_1e10/'];
path_out    = [parent_path, 'output/calvingMIP/Exp3-4/dx_2km/OceanVisc_1e10/'];

% Run in CECI cluster.
ctr.runmode  = 3;      % 1: graphics; 3: no graphics

%% Initial ice sheet creation
quarter   = true;
ctr.delta = 2e3;    % 5e3. Try run them at 2 km!!
ctr.imax  = 805;   % 322 (5 km), 805 (2 km), 1611 (1 km)
ctr.jmax  = 805;

Li    = (ctr.imax-1)*ctr.delta;
Lj    = (ctr.jmax-1)*ctr.delta;
[X,Y] = meshgrid(-Lj/2:ctr.delta:Lj/2,-Li/2:ctr.delta:Li/2);

R  = 800e3;   % 800e3
Bc = 900;     % 900
Bl = -2000;   % -2000
Ba = 1100;    % 1100
B  = BedGeom(X,Y,R,Bc,Bl,Ba);

ctr.m     = 3;
ctr.dt    = 1;
ctr.shelf = 1;
ctr.Asin  = zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up
ctr.Ao    = 2.9377e-18;

%Initial LSF mask
p       = nsidedpoly(1000, 'Center', [0 0], 'Radius', 755e3);
XX      = reshape(X, [numel(X),1]);
YY      = reshape(X, [numel(Y),1]);
IceMask = inpolygon(X,Y,p.Vertices(:,1),p.Vertices(:,2));

LSF             = zeros(ctr.imax,ctr.jmax);
LSF(IceMask==1) = 1;
LSF(IceMask==0) = -1;

% Create initial files.
out_1 = [path_in, 'ThuleLSF5'];
save(out_1, 'LSF');


%ctr.CF_Boundary = [path_in, 'CircThule.mat'];
ctr.LSFfile     = out_1;      % 'ThuleLSF5'

H  = zeros(ctr.imax,ctr.jmax)+10;
Mb = zeros(ctr.imax,ctr.jmax)+0.3;
Ts = zeros(ctr.imax,ctr.jmax)-5.0;

%---------------------------------------
% Cut out domain along symmetry axes (quarter).
if quarter == true

    ctr.mismip = 2;                 % Necessary for boundary conditions.
    ctr.imax   = (ctr.imax-1)/2+2;
    ctr.jmax   = (ctr.jmax-1)/2+2;
    ctr.Asin   = zeros(ctr.imax,ctr.jmax)+1e-7;

    B   = B(ctr.imax-2:end,ctr.jmax-2:end);
    H   = H(ctr.imax-2:end,ctr.jmax-2:end);
    Mb  = Mb(ctr.imax-2:end,ctr.jmax-2:end);
    Ts  = Ts(ctr.imax-2:end,ctr.jmax-2:end);
    LSF = LSF(ctr.imax-2:end,ctr.jmax-2:end);
end
%---------------------------------------

% Prepare initial file.
out_2 = [path_in, 'ThuleIn'];
save(out_2, 'B','H','Mb','Ts');


% Kori-ULB runs.
% 1. Initial spin up.
ctr.shelftune = 1;
ctr.SSA       = 1;         % ctr.SSA=1
ctr.dt        = 0.1;       % 0.2 (dx=2 km), 1, 2, 4. 
ctr.nsteps    = 60000;    % Jim: 10000; Daniel: 6000, 50000 (dt=0.2)
ctr.OceanVisc = 1e10;       % CalvingMIP: 8e9 (original), 7e9. Default 1e8.

out_3 = [path_out, 'Thule_quarter_Exp3_1'];
KoriModel(out_2, out_3, ctr);


% 2. Adjustment to imposition of Calving Front.
ctr.WV       = 0;
ctr.dt       = 0.1;     % 0.2, 1.0
ctr.calving  = 2;       % Direct, constant imposition of change in front positon.
ctr.LSFReset = 30;      % Jim: 50. Daniel: 30
ctr.nsteps   = 40000;   % Jim: 6000. Daniel: 4000 (dt=1). 20000

save(out_3,'LSF','-append');
out_4 = [path_out, 'Thule_quarter_Exp3_2'];
KoriModel(out_3, out_4, ctr); 


% 3. Impose zero rate of calving position change WV.
ctr.WV        = 0;
ctr.nsteps    = 1000;  % 100 (dt=1)

out_5 = [path_out, 'Thule_quarter_Exp3_3'];
%KoriModel(out_4, out_5, ctr); 


% 4. Calving MIP experiment 4 forcing.
ctr.calving   = 8;
ctr.nsteps    = 10000;   % 1000 (dt=1)
ctr.timeslice = 1;
ctr.LSFReset  = 50; %100, 250, 1000
ctr.snapshot  = 1000;
ctr.dt        = 0.1;  % Jim 0.05, Daniel: 1.0, 0.1
ctr.CR_AMP    = 750;

out_6 = [path_out, 'Thule_quarter_Exp4']
KoriModel(out_5, out_6, ctr); 


end

function [B]=BedGeom(x,y,R,Bc,Bl,Ba)
% param ters

rc=0;
%polarcoordinates
r=sqrt(x.*x+y.*y);
theta=atan2(y,x);
% B calculation
l=R-cos(2*theta).*R/2;
a=Bc-(Bc-Bl)*(r-rc).^2./(R-rc).^2;
B=Ba*cos(3*pi*r./l)+a ;
end