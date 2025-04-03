function grainsize = GrainSize(H,T,theta,EffStr,R,mu,D,M0,k0,gamma,p,c,n)

    % Daniel: grain size model. 
    % Based on Austin and Evans (2007), Ranganathan et al. (2021).

    % This model parametrizes migration recrystallization.
    % Glacial ice is a few tens of degrees below its melting 
    % temperature and thus deformation can warm ice to within a few degrees or 
    % less of its melting temperature. Migration of ice crystal becomes 
    % relevant in such cnoditions.
    % Grain size affect creep and vulnerability to ice fracture (damage).


    % define parameters for models
    H = 1000; % ice thickness, m
    %z = [0:10:H]; % m
    %N = 1000./10;
    n = 3; % flow law exponent
    p = 2; % grain-growth exponent
    c = 6; % for spherical grains (Behn et al in prep)
    Tm = 273; % melting temperature, in K
    Ts = Tm-25; % surface temperature, in K
    dT = Tm-Ts; % K
    rho = 917; % ice density, kg m^-3
    cp = 2050; % J kg^-1 K^-1
    K = 2.1; % W m^-1 K^-1
    Acons = 2.4e-24; % prefactor in Glen's flow law, Pa^-3 s^-1
    gamma = 0.065; % grain boundary energy, J/m^2 (Cuffey and Paterson 2010)
    k0 = 11.4266; %mm^p/s (Azuma et al 2012)
    k0 = k0./1000^p; % prefactor in the Arrhenius relation for grain growth, m^p/s
    R = 8.314; % J/mol K
    mu = 3e9; %Pa, shear modulus
    D = 0.05; %m, average grain size
    g = 9.8; % m s^-2
    M0 = 0.023; % intrinsic grain boundary mobility, m^2 s kg^-1
    theta = 0.99; % fraction of energy going into heating. 0.99
    
    % FIND GRAIN SIZE.
    % Activation energies.
    [Qg,Qc,Qm] = defineActivationEnergies(T);

    % Solve for strain rate.
    % input the Brinkmann number and Peclet number
    %Br = theta.*4;   % scale the Brinkmann number by fraction of energy going into heating
    %Pe = 2;          % Peclet number
    %strainrate = ((K.*dT.*Br)./(2.*H.^2.*theta)).^(n/(n+1)).*A.^(1./(n+1));

    % Factor.
    A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
    A = zeros(size(T));
    for i=1:length(T)
        A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
    end
    %tau = A.^(-1/n).*strainrate.^(1/n); % Pa
    tau = A.^(-1/n).*EffStr.^(1/n); % Pa

    % Compute grain size.
    grainsize = zeros(1,length(T));
    %grainsize = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T))) ./ ...
    %                (4.*mu.^2.*tau.*(1-theta).*strainrate)).^(1/(1+p));
    
    grainsize = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T))) ./ ...
                    (4.*mu.^2.*tau.*(1-theta).*EffStr)).^(1/(1+p));
    grainsize =  grainsize.*1e3;






end


