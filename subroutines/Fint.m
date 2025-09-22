function [F1x, F1y, F2, F2x, F2y] = Fint(ctr, eta_diva, H, zeta)
    
    % Preallocate arrays
    F1   = zeros([ctr.imax, ctr.jmax, ctr.kmax]);
    F2   = zeros([ctr.imax, ctr.jmax]);
    dz_H = zeros([ctr.imax, ctr.jmax, ctr.kmax]);
    z    = zeros([ctr.imax, ctr.jmax]);

    % Integrand order.
    n_2 = 2;


    % Differences in vertical dimension.
    % By definition: zeta(1)=0, zeta(ctr.kmax)=1;
    dz = zeros([ctr.kmax, 1]);
    dz(2:ctr.kmax) = diff(zeta);
    dz(1) = dz(2);


    zeta_m = 1.0 - zeta;


    % Vertical integration
    for k = 1:ctr.kmax

        % Redefine to adjust Kori def to integrals in Arthern et al. (2015).
        % Note that in Kori eta = eta * H.
        eta_diva(:,:,k) = eta_diva(:,:,k) ./ H;

        
        dz_H(:,:,k) = dz(k) * H;
        %z = z + zeta(k) * H;
        %H_minus_z = (H - z) ./ H;
        %value_1 = H_minus_z ./ eta_diva(:,:,k);
        %value_2 = H_minus_z.^ n_2 ./ eta_diva(:,:,k);

        value_1 = zeta_m(k) ./ eta_diva(:,:,k);
        value_2 = zeta_m(k).^ n_2 ./ eta_diva(:,:,k);

        
        %sum_1 = sum_1 + value_1;
        %F_1(:,:,k) = sum_1;
        F1(:,:,k) = F1(:,:,k) + value_1;
        F2        = F2 + value_2 * dz_H(k);
        
    end

    % Scale by dz_H.
    F1 = dz_H .* F1;

    % Stager since it is used to compute velocity and beta.
    F2x = 0.5 * ( F2 + circshift(F2, [0 -1]) );    % (i,j+1)
    F2y = 0.5 * ( F2 + circshift(F2, [-1 0]) );    % (i+1,j)

    F1x = 0.5 * ( F1 + circshift(F1, [0 -1 0]) );  % (i,j+1,k)
    F1y = 0.5 * ( F1 + circshift(F1, [-1 0 0]) );  % (i+1,j,k)

end