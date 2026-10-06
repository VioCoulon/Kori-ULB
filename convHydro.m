function convHydro(filePathIn, filePathOut, m1, m2, fl, swf, kappa, display)


%% Filenames
filePathTotoIn  = [filePathIn,  '_toto'];
filePathTotoOut = [filePathOut, '_toto'];

%% Parameters
nIter         = 1000;
coefRelax     = 0.5;
errorTol      = 1e-10;

%% Display
fprintf(['-------------------------------------',          '\n', ...
         'Converting initial state:',                      '\n', ...
         ' * inputFile  = ', filePathIn,                   '\n', ...
         ' * outputFile = ', filePathOut,                  '\n', ... 
         ' * initial state:',                              '\n', ...
         '    - power-law exponent: ', 'm1 = ', num2str(m1), '\n', ...
         '    - friction law:       ', fl,             '\n', ...
         '    - hydrology model:    ', 'NON',              '\n', ...
         ' * final state:',                                '\n', ...
         '    - power-law exponent: ', 'm2 = ', num2str(m2), '\n', ...
         '    - friction law:       ', fl                  '\n', ...
         '    - hydrology model:    ', 'NON',                '\n', ...
         '-------------------------------------',          '\n']);

%% Load everything from the NON file
load(filePathTotoIn);
par=KoriInputParams(ctr.m,ctr.basin);
ctr.mismip(any(ismember(fields(ctr),'mismip'))==0)=0;

%% Modify As from the NON state to get a swf state
ctr.m = m2;

switch(fl)
    case 'REGU'
        ctr.u0 = 300;
end

switch(swf)
    case 'NON'
        ctr.p            = 0;
        ctr.subwaterflow = 0;
    case 'HAB'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 0;
    case 'SWD'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 1;
    case 'TIL'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 2;
    case 'SWF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 3;
    case 'SCF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 4;
    case 'HARD'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 5;
    case 'SOFT'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 6;
    case 'HYB'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 7;
    case 'HARDINEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 8;
    case 'HARDEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 9;
    case 'SOFTINEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 10;
    case 'SOFTEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 11;
    case 'HYBINEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 12;
    case 'HYBEFF'
        ctr.p            = ctr.m;
        ctr.subwaterflow = 13;
end

% Target
AsfTarget = Asf;

% Set guess
AsGuess = Asf/par.AsScale;
As = AsGuess;
As0 = As;

% Particular cases
varAdd = {};
switch(swf)
    case 'SWD'
        varAdd = {'flw', 'Wd'};
    case 'SWF'
        varAdd = {'flw'};
    case 'TIL'
        Wtil = par.Wmax*ones(size(Asf));
        varAdd = {'flw', 'Wtil', 'Bmelt'};
    case 'SCF'
        varAdd = {'A', 'flw', 'ub', 'Bmelt'};
    case {'HARD', 'HARDINEFF','HARDEFF'}
        kappa = 0;
        varAdd = {'A', 'flw', 'ub', 'Bmelt'};
    case {'SOFT', 'SOFTINEFF', 'SOFTEFF'}
        kappa = 1;
        varAdd = {'A', 'flw', 'ub', 'Bmelt'};
    case {'HYB', 'HYBINEFF', 'HYBEFF'}
        varAdd = {'A', 'flw', 'ub', 'Bmelt', 'kappa'};
end

% Indices of values that must be changed
idx = (bMASK==0);

errorMean = Inf;
iter = 1;

ubSSA = ub;

while(errorMean >= errorTol)

    % Basal sliding
    [Asf,Asfx,Asfy,Asfd,Tbc,Neff,Wtil,r,expflw]=BasalSliding(ctr,par, ...
    As,Tb,H,B,MASK,Wd,Wtil,Bmelt,flw,bMASK,bMASKm,bMASKx,bMASKy);

    if(ctr.Tcalc >= 1)
        % SIA velocity
        [d,udx,udy,ud,ubx,uby,ub_new,uxsia,uysia,p,pxy]= ...
            SIAvelocity(ctr,par,A,Ad,Ax,Ay,Asfd,Asfx,Asfy,taud,G,Tb, ...
            H,Hm,Hmx,Hmy,gradm,gradmx,gradmy,gradxy,signx, ...
            signy,MASK,p,px,px,pxy);
        ub     = (1-coefRelax)*ub_new + coefRelax*ub;

        % Basal melt
        Bmelt_new = BasalMelting(ctr,par,G,taudxy,ub,H,tmp,dzm,MASK);
        Bmelt     = (1-coefRelax)*Bmelt_new + coefRelax*Bmelt;
    end

    % Subwater flux
    %[flw,Wd]=SubWaterFlux(ctr,par,H,HB,MASK,Bmelt);
    
    % Guess and target
    %par.effectHydroLimit = 5e-4;
    %effectHydro = max(par.effectHydroLimit, ((Neff/par.NeffScale).^(ctr.p))./expflw);
    effectHydro = 1.0;
    AsGuess =  Asf      ./(par.AsScale).*effectHydro;
    AsTarget = AsfTarget./(par.AsScale).*effectHydro;

    %if(strcmp(fl, 'REGU'))
    %    ussa     = vec2h(uxssa,uyssa);
    %    AsTarget = AsTarget ./ ((ussa + ctr.u0) ./ ctr.u0);
    %end

    ussa     = vec2h(uxssa,uyssa);

    if(strcmp(fl, 'WEERT-WEERT'))
        AsTarget = AsTarget.^(m2/m1) .* ussa.^(-m2*(m2-1)/(m1-1));
    end

    if(strcmp(fl, 'WEERT-REGU'))
        AsTarget = AsTarget.^(m2/m1) .* (ctr.u0 ./ (ussa+ctr.u0) ).^(-m2^2) .* ...
                                        ussa.^(-m2*(m2-1)/(m1-1));
    end

    % Compute mean error
    if(iter > 1)
        errorMean = mean(mean(abs((AsTarget(idx) - AsGuess(idx))./AsTarget(idx)*100)));
    end

    % Display
    if(display)
        disp(['Conversion of As (hydrology model) -- ', ....
            'iter ', num2str(iter), ': ', ...
            'mean relative error = ', num2str(errorMean), '%']);

        if(iter == 1)
            figure('Name', 'Conversion of As (hydrology model)');
        end
        semilogy(iter, errorMean, 'ob', 'linewidth', 1.5); hold on;
        grid on;
        xlabel('Iteration');
        ylabel('Mean relative error (%)');
        set(gca, 'Units', 'normalized', ...
            'FontUnits', 'points', ...
            'FontWeight','normal', ...
            'FontSize', 15, ...
            'FontName','CMU Sans Serif');
        drawnow;
    end
    
    % Iteration update
    if(iter >= nIter)
        warning(['The iterative scheme did not find a solution with the ', ...
                 'desired tolerance: ', ...
                 'error = ', num2str(errorMean), ' >= ', ...
                 'tol = ',   num2str(errorTol), '.']);
        break;
    else
        iter = iter + 1;
    end

    % Update
    As(idx) = As(idx) + (1-coefRelax) * (AsTarget(idx) - AsGuess(idx));
end
ub = ubSSA;

%% Save summary file
varSumm = fieldnames(load(filePathIn));
save(filePathOut, varSumm{1});
for i = 2:length(varSumm)
    var = varSumm{i};
    save(filePathOut, var, '-append');
end
for i = 1:length(varAdd)
    var = varAdd{i};
    save(filePathOut, var, '-append');
end

%% Save complete file
varFull = fieldnames(load(filePathTotoIn));
save(filePathTotoOut, varFull{1});
for i = 2:length(varFull)
    var = varFull{i};
    save(filePathTotoOut, var, '-append');
end
for i = 1:length(varAdd)
    var = varAdd{i};
    save(filePathTotoOut, var, '-append');
end


%% Plots.
% Original with m1.
myfig(log10(As0),-9,-1)


% New with m2.
myfig(log10(As),-9,-1)


% Difference.
myfig(log10(abs(As-As0)),-9,-1)


end