function noise = stochastic_noise(ctr, sigma, tau)
    % Stochastic noise function.
    % Parameters:
    % dt   - Time step.
    % sigma - Standard deviation of noise.
    % tau  - Persistence time.
    % N    - Number of points.

    % For odd number of time steps, make it even by increasing 1.
    if mod(ctr.nsteps,2) == 0
        N = ctr.steps;
    else
        N = ctr.nsteps + 1;
    end

    % Frequency resolution and maximum frequency.
    df = 1.0 / (ctr.dt * N);
    f0 = 0.5 / ctr.dt;

    % Frequency array.
    f1 = 0:df:(f0 - df);

    % Auto-correlation at a lag of dt.
    r = 1.0 - (ctr.dt / tau);

    % Scales total variance.
    P_0 = 1.0;

    % Analytical power spectra from a given auto-correlation r.
    P = sqrt(P_0 ./ (1.0 + r^2 - 2.0 * r * cos(2.0 * pi * ctr.dt * f1)));

    % Seed random number generator.
    rng(1);

    % Create array with random phase.
    phase_all = 1i * 2.0 * pi * rand(1, N);

    % Concatenate power spectra with different persistence time tau.
    P_r2 = [P, fliplr( P(1:floor(0.5*N)) )];

    %fprintf('\n P_r2 = ', P_r2);
    %length(P_r2)
    %fprintf('\n phase_all = ', phase_all);
    %length(phase_all)

    % Include identical random phase.
    P_rand = P_r2 .* exp(phase_all);

    % Transform back to time space to save stochastic signal.
    noise = normalised_ifft(P_rand, sigma);
end

function output = normalised_ifft(P_rand, sigma)
    
    % Perform inverse FFT and normalize the signal.
    ifft_result = ifft(P_rand, 'symmetric');
    
    % Normalize the result by scaling with the desired standard deviation.
    %output = sigma * ifft_result / std(ifft_result);
    output = sigma * ( ifft_result - mean(ifft_result) ) / std(ifft_result);
end