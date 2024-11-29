% Basic tests for the Kori Model 
% tests to be carried out to test code when changes occured



function BasicTests(n)

    close all;
    
    scenario_txt={'All tests','Mass conservation','EISMINT', ...
        'MISMIP circular','Circular quarter', ...
        'MISMIP+','Antarctica init SIA','Antarctica Init Hybrid', ...
        'Antarctica Run','Aletsch glacier','Basin init SIA', ...
        'Basin init Hybrid','Basin Run'};
    if nargin<1
        fprintf('Experiment:\n');
        for i=1:length(scenario_txt)
            fprintf('  (%d) %s\n',i-1,char(scenario_txt(i)));
        end
        n=input('Choose scenario: ');
    end
    switch n
        case 0
            mascon;
            EismintTest;
            MismipTest;
            Circular;
            MismipPlus;
            Antarctica(1);
            Antarctica(2);
            Antarctica(3);
            Aletsch;
            RunASE(1);
            RunASE(2);
            RunASE(3);
        case 1
            mascon;
        case 2
            EismintTest;
        case 3
            MismipTest;
        case 4
            Circular;
        case 5
            MismipPlus;
        case 6
            Antarctica(1);
        case 7
            Antarctica(2);
        case 8
            Antarctica(3);
        case 9
            Aletsch;
        case 10
            RunASE(1);
        case 11
            RunASE(2);
        case 12
            RunASE(3);
        otherwise
            disp('wrong value')
    end

end

function mascon

    ctr.imax=61;
    ctr.jmax=61;
    ctr.delta=25.e3;
    ctr.nsteps=501;
    ctr.dt=10.;
    ctr.Asin=5e-9; %3e-9
    ctr.SSA=2;
    ctr.ItSolv=0; % makes a difference
    ctr.upstream=0; % makes a big difference: 0 more precise

    H=zeros(ctr.imax,ctr.jmax);
    start=(ctr.imax-1)/2+1-10;
    stop=(ctr.imax-1)/2+1+10;
    H(start:stop,start:stop)=1000;
    save('MASCON','H');

    KoriModel('MASCON','mascon1',ctr);


end

function EismintTest

    % TEST 1: EISMINT
        
    ctr.imax=31; %31
    ctr.jmax=31; %31
    ctr.kmax=21;
    
    ctr.delta=50.e3;
    ctr.nsteps=5001; % 2001
    ctr.snapshot=11;
    ctr.dt=5; % 50
    ctr.Tcalc=1;
    ctr.plotH=1;
    ctr.Asin=0;
    ctr.MbConst=0.3;
    
    Li=(ctr.imax-1)*ctr.delta/1e3;
    Lj=(ctr.jmax-1)*ctr.delta/1e3;
    [X,Y] = meshgrid(0:ctr.delta/1e3:Li,0:ctr.delta/1e3:Lj);
    dist=max(abs(X-Lj/2.),abs(Y-Li/2.));

    % EISMINT fixed margin
    Ts=239.+8e-8*dist.^3.-273.15;
    save('EismintIn','Ts');
    KoriModel('EismintIn','Eismint1',ctr);
    PlotFiguresEISMINT(1,'Eismint1');

    % EISMINT moving margin
    ctr.MbType=4;
    ctr.TsType=2;
    KoriModel('EismintIn','Eismint2',ctr);
    PlotFiguresEISMINT(2,'Eismint2');
    
end


function PlotFiguresEISMINT(type,output)

    filename=[output,'_toto'];
    load(filename);
    iplot=(ctr.imax+1)/2;

    if type==1
        figure;
        Temp=Tb+par.pmp*H;
        plot(x(iplot:end)-750,Temp(iplot,iplot:end),'o-'); hold on;
        dataX=[0 50 100 150 200 250 300 350 400 450 500 550 600 650 700];
        dataT=[-9 -8.5 -8 -7 -6 -5 -4 -2.5 -0.5 0 0 0 0 0 0];
        plot(dataX,dataT,'xr-');
        xlabel('Distance from center (km)');
        ylabel('Homologous basal temperature (�C)');
        grid on;
        xlim([0 700]);
        
        fprintf('EISMINT FM\nH divide = %f\n',max(max(H))); % ice thickness in center
        % mass flux at midpoint
        midp=iplot/2;
        dmean=(d(iplot,midp)+d(iplot,midp-1)+d(iplot-1,midp)+d(iplot-1,midp-1))/4;
        slope=(H(iplot,midp+1)-H(iplot,midp-1))/(2*ctr.delta);
        fprintf('Flux midpoint = %f\n',dmean*slope);
        % divide basal temperature (homologuous)
        fprintf('T divide = %f\n\n',Temp(16,16));
    else
        figure;
        Temp=Tb+par.pmp*H;
        Temp(H<=5)=0;
        plot(x(iplot:end)-750,Temp(iplot,iplot:end),'o-'); hold on;
        dataX=[0 50 100 150 200 250 300 350 400 450 500 550 600 650 700];
        dataT=[-13.4 -12.3 -11 -9.6 -8 -6 -3.7 -1 0 0 0 0 0 0 0];
        plot(dataX,dataT,'xr-');
        xlabel('Distance from center (km)');
        ylabel('Homologous basal temperature (�C)');
        grid on;
        xlim([0 700]);
        
        fprintf('EISMINT MM\nH divide = %f\n',max(max(H))); % ice thickness in center
        % mass flux at midpoint
        midp=iplot/2;
        dmean=(d(iplot,midp)+d(iplot,midp-1)+d(iplot-1,midp)+d(iplot-1,midp-1))/4;
        slope=(H(iplot,midp+1)-H(iplot,midp-1))/(2*ctr.delta);
        fprintf('Flux midpoint = %f\n',dmean*slope);
        % divide basal temperature (homologuous)
        fprintf('T divide = %f\n\n',Temp(16,16));
    end

end


function MismipTest
    % TEST 2: Circular MISMIP ice sheet
    % Test on Schoof and symmetry of ice sheet
    
    ctr.schoof=1;
    ctr.imax=67;
    ctr.jmax=67;
    ctr.delta=50.e3;
    ctr.m=2;
    ctr.nsteps=1001; %2001, 10001
    ctr.dt=10;       % 5, 1
    ctr.SSA=2;
    ctr.shelf=1;
    ctr.shelftune=ones(ctr.imax,ctr.jmax);
    
    Li=(ctr.imax-1)*ctr.delta/1e3;
    Lj=(ctr.jmax-1)*ctr.delta/1e3;
    [X,Y] = meshgrid(0:ctr.delta/1e3:Li,0:ctr.delta/1e3:Lj);
    dist=sqrt((X-Lj/2.).^2.+(Y-Li/2.).^2.);
    B=720-778.5*dist/750.;
    H=zeros(ctr.imax,ctr.jmax);
    Mb=zeros(ctr.imax,ctr.jmax)+0.3;
    Ts=zeros(ctr.imax,ctr.jmax)-10.;
    ctr.Asin=zeros(ctr.imax,ctr.jmax)+3.0e-9;
    
    %Initial LSF mask (R2016b)
    LSF=ones(ctr.imax,ctr.jmax);
    LSF(dist>1500)=-1;
    H(LSF<0)=0;
    ctr.calving=2;
    ctr.WV=0;

    save('MismipIn','B','H','Mb','Ts','LSF');

    % Daniel.
    ctr.timeslice=1;
    ctr.snapshot=25; % 800 output every 10 years.

    KoriModel('MismipIn','mismip1',ctr);
     
end


function Circular

    % Initial ice sheet creation

    ctr.delta=10e3;
    ctr.imax=161; 
    ctr.jmax=161;

    Li=(ctr.imax-1)*ctr.delta;
    Lj=(ctr.jmax-1)*ctr.delta;
    [X,Y]=meshgrid(-Lj/2:ctr.delta:Lj/2,-Li/2:ctr.delta:Li/2);

    R=800e3 ;%800e3
    Bc=900;  %900
    Bl=-2000; %-2000
    B=BedGeom(X,Y,R,Bc,Bl);

    %Initial LSF mask (R2016b)
    numSides = 1000;
    center = [0 0];
    radius = 755e3;
    theta = linspace(0, 2*pi, numSides+1);
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    IceMask = inpolygon(X, Y, x, y);
    LSF=zeros(ctr.imax,ctr.jmax);
    LSF(IceMask==1)=1;
    LSF(IceMask==0)=-1;

    %---------------------------------------
    % Cut out domain along symmetry axes
    ctr.imax=(ctr.imax-1)/2+2;
    ctr.jmax=(ctr.jmax-1)/2+2;
    B=B(ctr.imax-2:end,ctr.jmax-2:end);
    LSF=LSF(ctr.imax-2:end,ctr.jmax-2:end);
    ctr.mismip=2;
    %---------------------------------------

    ctr.m=3;
    ctr.dt=1;
    ctr.shelf=1;
    ctr.shelftune=1;
    ctr.SSA=1;

    ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-7; % Same as Hilmars set up
    ctr.Ao=2.9377e-18;

    H=zeros(ctr.imax,ctr.jmax)+10;
    Mb=zeros(ctr.imax,ctr.jmax)+0.3;
    Ts=zeros(ctr.imax,ctr.jmax)-5.0;
    H(LSF<0)=0.1;
    MASKo=ones(ctr.imax,ctr.jmax);
    MASKo(LSF<0)=0;
    save('CircularIn','B','H','Mb','Ts','LSF','MASKo');

    %Iinitial spin up
    ctr.nsteps=1001;
    ctr.calving=2;
    ctr.WV=0;
    % Daniel.
    ctr.timeslice=1;
    ctr.snapshot=25; % 800 output every 10 years.
    KoriModel('CircularIn','Circular1',ctr); 

end


function [B]=BedGeom(x,y,R,Bc,Bl)
    % parameters

    rc=0;
    %polarcoordinates
    r=sqrt(x.*x+y.*y);
    theta=atan2(y,x);
    % B calculation
    l=R-cos(2*theta).*R/2;
    a=Bc-(Bc-Bl)*(r-rc).^2./(R-rc).^2;
    B=a ;
end



function MismipPlus

% MISMIP+ experiment with Kori-ULB
% Details in Cornford et al (2020)

% Initial ice sheet creation

ctr.delta=2e3;
ctr.imax=23; % need number of cells + 2
ctr.jmax=352; % need number of cells + 1
ctr.m=3;
ctr.dt=1;
ctr.nsteps=5001;
ctr.mismip=1; % indicates that BCs need to be applied
ctr.SSA=1;
ctr.shelf=1;
ctr.Ao=4e-17; % 4.0e-17
ctr.shelftune=1;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-6; % 1e-5

Li=(ctr.imax-2)*ctr.delta;
Lj=(ctr.jmax-2)*ctr.delta;
[X,Y]=meshgrid(-ctr.delta:ctr.delta:Lj,-ctr.delta:ctr.delta:Li);
B=Bx(X)+By(Y,4.0e3,5.0e2,24.0e3);
B=max(B,-720.0);
% periodic BCs
B(1,:)=B(3,:);
B(ctr.imax,:)=B(ctr.imax-2,:);
H=zeros(ctr.imax,ctr.jmax)+10;
Mb=zeros(ctr.imax,ctr.jmax)+0.3;
Ts=zeros(ctr.imax,ctr.jmax)-5.0;
save('MismipPlusIn','B','H','Mb','Ts');

KoriModel('MismipPlusIn','MismipPlus1',ctr);

% Ice1rr experiment for 200 year with melting

ctr.nsteps=201;
ctr.meltfunc=10;
ctr.meltfac=1;
ctr.BetaIter=ctr.nsteps;
KoriModel('MismipPlus1','MismipPlus2',ctr);

% Ice1ra experiment starting from Ice1r without melt

ctr.nsteps=101;
ctr.meltfunc=0;
KoriModel('MismipPlus2','MismipPlus3',ctr);

end


function [Bxx]=Bx(x)
    B0=-150;
    B2=-728.8;
    B4=343.91;
    B6=-50.57;
    xx=x./300.0e3;
    xx2=xx.^2;
    xx4=xx2.^2;
    xx6=xx4.*xx2;
    Bxx=B0+B2*xx2+B4*xx4+B6*xx6;
end


function [Byy]=By(y,fc,dc,wc)
    Byy=dc./(1.0+exp(-2.0*(y-wc)/fc)) + dc./(1.0+exp(2.0*(y+wc)./fc));
end


function Antarctica(n)

    % whole continent Antarctica
    % inversion 1 and 2 + run
    
    ctr.inverse=1;
    ctr.imax=225;
    ctr.jmax=225;
    ctr.kmax=21;
    ctr.delta=25.e3;
    ctr.nsteps=1001;
    ctr.dt=10;
    ctr.Tcalc=2;
    ctr.Tinit=1; % SET 1 FOR INITIALIZATION !

    ctr.Asin=3e-9; %3e-9
    ctr.Ao=5.0e-17; % 5.0e-17
    ctr.m=3;
    ctr.calving=2;
    ctr.WV=0;

    if n==1
        KoriModel('Bedmachine25km','INIT25a',ctr);
    end

    ctr.inverse=2;
    ctr.Tinit=0;
    ctr.SSA=2;
    ctr.shelf=1;
    ctr.schoof=1; % 1
    ctr.dt=0.2; % 0.1
    ctr.nsteps=201;
    ctr.Tinv=20;

    if n==2
        KoriModel('INIT25a','INIT25b',ctr);
    end

    ctr.inverse=0;
    ctr.PDDcalc=1;
    ctr.meltfunc=3;
    ctr.calving=4;
    fc.DeltaT=zeros(ctr.nsteps,1)+10;
    
    if n==3
        KoriModel('INIT25b','Run25a',ctr,fc);
    end
    
end
  

function Aletsch

% test on modelling Aletsch glacier
% assembling input files


ctr.imax=157;
ctr.jmax=113;
ctr.delta=150;
ctr.nsteps=701;
ctr.dt=0.5;
ctr.plotH=2;

ctr.MbType=6; % 6 McCall Mb
ctr.TsType=3; % 3 McCall Ts
ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-9;
ctr.ELA=3000;
ctr.mbgrad=0.008;
ctr.mbgrad1=0.003;

fc.DeltaT=zeros(ctr.nsteps,1);
ctr.Tcalc=1;
KoriModel('Aletsch150m','Aletsch1',ctr,fc);


end


function RunASE(n)

% test on modelling PIG and Thwaites at 3km
% assembling input files
% Kori-ULB

% First SIA inversion (As)
ctr.imax=318;
ctr.jmax=270;
ctr.kmax=21;
ctr.delta=3.e3;
ctr.nsteps=601;
ctr.dt=5;
ctr.Ao=5.0e-17; % 5.0e-17
ctr.m=3;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.calving=2;
ctr.WV=0;
ctr.Tinit=1;
ctr.Tcalc=2;
ctr.basin=1;
ctr.inverse=1;

if n==1
    KoriModel('ASE3km','ASEint1',ctr);
end


% Second SSA inversion (As, MeltInv)

ctr.inverse=2;
ctr.meltfunc=1;
ctr.GroundedMelt=1; % better for basins!
ctr.shelf=1;
ctr.SSA=2;
ctr.Tinit=0; % SET 1 FOR INITIALIZATION !
ctr.Tinv=10;
ctr.nsteps=251;
ctr.dt=0.2;

if n==2
    KoriModel('ASEint1','ASEint2',ctr);
end

% Forcing run
ctr.inverse=0;
ctr.calving=4;
ctr.meltfunc=3;
ctr.gammaT=1e-3;

if n==3
    KoriModel('ASEint2','ASErun1',ctr); % start from optimized run
end


end

