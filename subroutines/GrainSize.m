function d_grain = GrainSize(H,T,EffStr,zeta,ctr,par)

    % Daniel: grain size model. 
    % Based on Austin and Evans (2007), Ranganathan et al. (2021).

    % This model parametrizes migration recrystallization.
    % Glacial ice is a few tens of degrees below its melting 
    % temperature and thus deformation can warm ice to within a few degrees or 
    % less of its melting temperature. Migration of ice crystal becomes 
    % relevant in such cnoditions.
    % Grain size affect creep and vulnerability to ice fracture (damage).


    % define parameters for models
    %H = 1000; % ice thickness, m
    %z = [0:10:H]; % m
    %N = 1000./10;
    n = 3; % flow law exponent
    p = 2; % grain-growth exponent
    c = 6; % for spherical grains (Behn et al in prep)
    Tm = 273; % melting temperature, in K
    %Ts = Tm-25; % surface temperature, in K
    %dT = Tm-Ts; % K
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

    %min(T(:))
    %max(T(:))

    %min(EffStr(:))
    %max(EffStr(:))
    

    % Solve for strain rate.
    % input the Brinkmann number and Peclet number
    %Br = theta.*4;   % scale the Brinkmann number by fraction of energy going into heating
    %Pe = 2;          % Peclet number
    %strainrate = ((K.*dT.*Br)./(2.*H.^2.*theta)).^(n/(n+1)).*A.^(1./(n+1));

    % Factor.
    %A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
    %A = zeros(size(T));
    %for i=1:length(T)
    %    A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
    %end



    % THERE IS A PROBLEM WHEN USING THE TEMPERATURE INPUT FROM KORI!!!!!
    % CONSTANT TEMPERATURES TEST.
    %T(:) = 265.0;
    %T(T<260) = 260.0;

    % Activation energies.
    [Qg,Qc,Qm] = ActivationEnergies(T);


    
    % Ranganathan et al., (2021).
    %A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
    %A = A0*exp(-(Qc/R).*((1/T)-(1/263)));

    % Greve and Blatter (2009)????
    A = 1.916e3 * exp(-139000/(R*T));

    % Test units.
    %A = A * par.secperyear;
    
    % Ranganathan.
    %tau = A.^(-1/n).*strainrate.^(1/n); % Pa

    % In their model, strain rates are expressed in s^-1.
    %EffStr =  EffStr * par.secperyear;
    %EffStr =  EffStr / par.secperyear;

    % Strain rate definition as in Moreno-Parada et al. (2024).
    %strainrate = 2.0 * A .^ (-1.0/n) .* EffStr .^((n+1)/n);
    
    % This seem to work? Double check!!!!!!!!!!!!!!!!!!!!!!!!!
    EffStr = EffStr.^(1.0/n);
    %EffStr = EffStr.^((1-par.n)/(2*par.n));       % Kori definition.
    %EffStr = EffStr.^(1.0/(2*n));
    EffStr_3D = repmat(EffStr, [1, 1, ctr.kmax]);

    % Test 1.
    tau = A.^(-1.0/n) .* EffStr_3D; % Pa


    % Compute grain size.
    %grainsize = zeros(1,length(T));

    % Ranganathan.
    %grainsize = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T))) ./ ...
    %                (4.*mu.^2.*tau.*(1-theta).*strainrate)).^(1/(1+p));
    
    % Daniel.
    d_grain_3D = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma + tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T))) ./ ...
                    (4.*mu.^2.*tau.*(1-theta).*EffStr_3D)).^(1/(1+p));


    %min(d_grain_3D(:))
    %max(d_grain_3D(:))
    
    % They express it in milimetres, but we need metres.
    %d_grain_3D =  d_grain_3D.*1e3;

    % Grain size must be transformed into a 2D field.
    
    % Reshapes and replicates the vector dz into a 3D array.
    dz = zeros([ctr.kmax, 1]);
    dz(2:ctr.kmax) = diff(zeta);
    dz(1) = dz(2);
    zl = repmat(reshape(dz, 1, 1, ctr.kmax), [ctr.imax, ctr.jmax, 1]);
    
    %size(dz)
    %size(zl)
    %size(d_grain_3D)

    % Vertical average considering unevenly-spaced vertical grid in Kori.
    %weighted_sum = sum(d_grain_3D .* zl, 3);    % Element-wise multiply and sum along k.
    %sum_weights  = sum(zl, 3);                  % Sum of weights along k.
    %d_grain = weighted_sum ./ sum_weights;      % Compute weighted mean to obtain vertically averaged 2D field.

    % Regular average.
    %d_grain = mean(d_grain_3D(:,:,1:5), 3);


    % Take a certain level. k = 11 is the base.
    d_grain = d_grain_3D(:,:,11);

    d_grain(d_grain<1.0e-5) = 1.0e-5;


end


