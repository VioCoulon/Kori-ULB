function convPowerLaw(filePathIn, filePathOut, mIn, mOut)

%% Path
%if(isInFolder('init'))
%    addpath(genpath('../'));
%end

%% Filenames
%filePathTotoIn  = [filePathIn,  '_toto'];
%filePathTotoOut = [filePathOut, '_toto'];


filePathTotoIn  = [remMatExt(filePathIn),  '_toto'];
filePathTotoOut = [remMatExt(filePathOut), '_toto'];

filePathIn      = remMatExt(filePathIn);
filePathTotoIn  = remMatExt(filePathTotoIn);
filePathOut     = remMatExt(filePathOut);
filePathTotoOut = remMatExt(filePathTotoOut);

%% Checks
% Check if input files exist
if(not(exist(addMatExt(filePathIn), 'file') == 2))
    error(['The file <', filePathIn,     '> does not exist.']);
end
if(not(exist(addMatExt(filePathTotoIn), 'file') == 2))
    error(['The file <', filePathTotoIn, '> does not exist.']);
end

%% Parameters
AsMax = 1e+03;
AsMin = 1e-20;

%% Display
fprintf(['-------------------------------------',          '\n', ...
         'Converting power-law initial state:',            '\n', ...
         ' * inputFile  = ', filePathIn,                   '\n', ...
         ' * outputFile = ', filePathOut,                  '\n', ...
         ' * initial state:',                              '\n', ...
         '    - source:      m = ', num2str(mIn),          '\n', ...
         '    - destination: m = ', num2str(mOut),         '\n', ...
         '-------------------------------------',          '\n']);


%% Load everything from the m = mIn file
data      = load(filePathIn);
data_toto = load(filePathTotoIn);

% Assign data.
u   = data_toto.u;
As  = data_toto.As;

% Define u0 in Coulomb law if starting from Weertman.
u0  = 300.0;


As0 = As;

%% Modify As from the m = mIn state to get a m = mOut state
AsScaleIn  = 1e5^(2-mIn);
AsScaleOut = 1e5^(2-mOut);

% From Weert to Weert with different exponents.
%As = (1/AsScaleOut) * (As * AsScaleIn).^(mOut/mIn) .* u.^(1 - mOut/mIn);

% From Weertman to Coulomb law.
As = (1/AsScaleOut) * (As * AsScaleIn).^(mOut/mIn) .* u.^(1 - mOut/mIn) .* (u0./(u+u0));

As = min(As, AsMax);
As = max(As, AsMin);

Asf = AsScaleOut * As;

%% Save summary file
%varSumm = fieldnames(load(filePathIn));
%save(filePathOut, varSumm{1});
%for i = 2:length(varSumm)
%    var = varSumm{i};
%    save(filePathOut, var, '-append');
%end

data.As          = As;
%data.Asf         = Asf;
%data.ctr.m       = mOut;
%data.par.AsScale = AsScaleOut;

save(filePathOut, '-struct', 'data');



%% Save complete file
data_toto.ctr.m       = mOut;
data_toto.par.AsScale = AsScaleOut;

%varFull = fieldnames(load(filePathTotoIn));
%save(filePathTotoOut, varFull{1});
%for i = 2:length(varFull)
%    var = varFull{i};
%    save(filePathTotoOut, var, '-append');
%end

data_toto.As          = As;
data_toto.Asf         = Asf;
data_toto.ctr.m       = mOut;
data_toto.par.AsScale = AsScaleOut;



save(filePathTotoOut, '-struct', 'data_toto');


function strOut = addMatExt(strIn)

    ext = '.mat';

    if(length(strIn) > 4)
        if(strcmp(strIn(end-3:end), '.mat'))
            ext = [];
        end
    end

    strOut = [strIn, ext];
end

function strOut = remMatExt(strIn)

    len = length(strIn);

    if(len > 4)
        if(strcmp(strIn(end-3:end), '.mat'))
            len = len-4;
        end
    end

    strOut = strIn(1:len);
end


%% Plots.
% Original with m1.
myfig(log10(As0),-10,-1)
title('Initial As (m1)');

% New with m2.
myfig(log10(As),-10,-1)
title('Final As (m2)');

% Difference.
%myfig(log10(abs(As-As0)),-10,-1)


end