

% Matlab script to edit Kori run files varying certain parameter.
% Daniel Moreno Parada - March 2025.
% daniel.moreno.parada@ulb.be


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORWARD RUN (CONSTANT FORCING).
ctr.imax     = 472;
ctr.jmax     = 400;
ctr.delta    = 2.e3;
ctr.Ao       = 5.0e-17; % 5.0e-17
ctr.m        = 3;
ctr.Asin     = zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.basin    = 1;    % Run the model for a specific basin.

ctr.inverse    = 0;
ctr.shelf      = 1;        % Ice shelves are considered.
ctr.SSA        = 2;
ctr.calving    = 4;          % 5, to apply LSF function.
ctr.dt         = 0.05;        % 0.02, 0.1, 0.05, 0.025
ctr.nsteps     = 15001;       % 20001
ctr.meltfunc   = 3;          % PICO
ctr.gammaT     = 10.0e-5;     % 2.5e4, 1.0e-4, 0.25e-4. Best fit in PICO: 2e-5.
ctr.meltfac    = 1;          % Factor multiplying sub-shelf melt.
ctr.LimitFront = 1;          % Calving front limit by initial location.

ctr.timeslice = 1;
ctr.snapshot  = 50;          % 50, 1500 (dt=0.1 yr).

ctr.stochastic = 0;
ctr.sigma_To   = 0.5;
ctr.sigma_Mb   = 0.3;
ctr.tau_Mb     = 1.0;
ctr.tau_To     = 10.0;
ctr.seed       = 100;



% Define different values for ctr.gammaT
gammaT_values = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5];

script_names = cell(length(gammaT_values), 1);


% Experiment and variables to be changed.
%name_1 = 'sigma_oce';
%name_2 = 'seed';

% Deterministic.
name_1 = 'deter';
name_2 = 'HR';


sigma_now = ctr.sigma_To*100;
if sigma_now < 100
    num = sprintf('0%.0f', sigma_now);
    num
else
    num = sprintf('%.0f', sigma_now);
    num
end


% Stochastic.
exp_1 = sprintf('%s%s', name_1, num);
%exp_2 = sprintf('%s%d', name_2, ctr.seed);   % Stoch exps.

% Deterministic.
exp_1 = name_1;
exp_2 = 'HR';                              % Deter exps.

var   = 'gamma';

% Full paths.
% Nic5.
exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/nic5/';
parent_path    = '/scratch/ulb/glaciol/dmoreno/Kori-ULB/';
rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/stoch/';                             % Stoch.
rel_out        = 'output/thwaites/deter/';                              % Deter. 

% Local.
%parent_path = '/home/daniel/models/Kori-ULB/';
%rel_in      = 'output/Thwaites/Initialization/eta1e7/ground_melt_0/';
%rel_out     = 'output/Thwaites/Stochastic/test/';

path_scripts  = [exe_path_local, exp_1, '/'];
path_new_file = [path_scripts, exp_2, '/'];



% Loop to generate filenames
for i = 1:length(gammaT_values)
    value_now = 1e5*gammaT_values(i);
    if i<10
        num = sprintf('0%.0fe5', value_now);
    else
        num = sprintf('%.0fe5', value_now);
    end

    script_names{i} = [exp_2, '_', var, num];
    
end

% Read the original file
original_file = 'RunASE_nic5.m'; 
file_contents = fileread(original_file);


% Loop through each gammaT value and save a new file
for i = 1:length(gammaT_values)

    % Value gamma.
    new_value = sprintf('ctr.gammaT = %.5e;', gammaT_values(i));
    
    % Create the folder if it does not exist
    folder_name = [path_new_file, script_names{i}]
    
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end

    % Identical param file name for consistency when reading.
    full_mat = [folder_name, '/params.mat']

    % Update control values.
    ctr.gammaT = gammaT_values(i);

    % We update the name of the file in accordance with the value of gammaT.
    % Avoid sign "-" in the file name as it does not compile.
    % Appropirate numering to ensure order when listing in Linux.
    value_now = 1e5*ctr.gammaT;
    if value_now < 10
        num = sprintf('0%.0fe5', value_now);
    else
        num = sprintf('%.0fe5', value_now);
    end

    % Input and output file names.
    name_1 = 'INIT_SSA_3'; % ASEfor_m2_80, INIT_NON_3, INIT_NON_LONG_10K
    name_2 = [exp_2, '_', var, num];

    path_in     = [parent_path, rel_in];
    %path_out    = [parent_path, rel_out, exp_1, '/', exp_2, '/', name_2, '/'];
    path_out    = [parent_path, rel_out, exp_1, '/', name_2, '/'];

    input       = [path_in, name_1];
    output      = [path_out, name_2];

    % Save all variables to a .mat file
    save(full_mat, 'ctr', 'input', 'output');
    
    fprintf('Saved: %s\n', full_mat);
    fprintf('Input: %s\n', input);
    fprintf('Output: %s\n', output);
end