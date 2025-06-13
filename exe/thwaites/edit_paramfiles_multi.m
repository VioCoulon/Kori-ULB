

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

ctr.runmode    = 3;  % 1: graphics; 3: no graphics
ctr.inverse    = 0;
ctr.shelf      = 1;        % Ice shelves are considered.
ctr.SSA        = 2;
ctr.calving    = 4;          % 5, to apply LSF function.
ctr.Tcalc      = 2;            % Calculate temperature field and thermomechanical coupling, , i.e. A = f (T)
ctr.dt         = 0.05;        % 0.02, 0.1, 0.05, 0.025
ctr.nsteps     = 20001;       % 15001, 20001
ctr.meltfunc   = 3;          % PICO
ctr.gammaT     = 2.0e-5;     % 2.5e4, 1.0e-4, 0.25e-4. Best fit in PICO: 2e-5.
ctr.meltfac    = 1;          % Factor multiplying sub-shelf melt.
ctr.LimitFront = 1;          % Calving front limit by initial location.

ctr.timeslice = 1;
ctr.snapshot  = 50;          % Normal: 50. HR: 1500 (dt=0.1 yr).

ctr.stochastic = 1;      % Deter: 0. Stoch: 1.
ctr.sigma_To   = 4.0;    % 0.1, 0.25, 0.5, 1.0, 2.0, 4.0
ctr.sigma_Mb   = 0.3;
ctr.tau_Mb     = 1.0;
ctr.tau_To     = 20.0;   % 1, 5, 10, 20, 70 yr.
ctr.seed       = 100;


% Define different values for ctr.gammaT.
values_1 = [0.125, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
%values_1 = [0.5, 1.5, 2.5, 3.0, 3.5, 4.5, 5];

%values_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
values_2 = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29];
%values_2 = [0];

% Create empty dictionary to populate it with exp names.
script_names = cell(length(values_1), length(values_2));

% Experiment name.
exp_1 = 'sigma_oce400';
%exp_1 = 'deter';

% Variables.
var_1   = 'meltfac';      % gamma
var_2   = 'seed';     % snapshot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%25%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL PATHS.

% Nic5.
%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/nic5/stoch/tau_To_70/';
%parent_path    = '/scratch/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/stoch/tau_To_70/';                             % Stoch.

%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/nic5/';
%parent_path    = '/scratch/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/';  


% Lyra.
%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lyra/stoch/tau_To_01/';
%parent_path    = '/globalsc/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/stoch/tau_To_01/';                             % Stoch.


% lemaitre4.
exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lemaitre4/stoch/extra_runs/tau_To_70/';
parent_path    = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/';
rel_in         = 'ice_data/eta1e7/ground_melt_0/';
rel_out        = 'output/thwaites/stoch/extra_runs/tau_To_70/';                             % Stoch.

%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lemaitre4/';
%parent_path    = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/deter/';                             % Stoch.


% Local.
%parent_path = '/home/daniel/models/Kori-ULB/';
%rel_in      = 'output/Thwaites/Initialization/eta1e7/ground_melt_0/';
%rel_out     = 'output/Thwaites/Stochastic/test/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



path_scripts  = [exe_path_local, exp_1, '/'];



% Loop to generate filenames
for i = 1:length(values_1)
    for j = 1:length(values_2)

        % Gamma
        % value_now = 1e5*values_1(i);
        %if value_now<10
        %    num_1 = sprintf('0%.0fe5', value_now);
        %else
        %    num_1 = sprintf('%.0fe5', value_now);
        %end

        % Melt fac.
        value_now = 1e3*values_1(i);

        if value_now<1e3
            num_1 = sprintf('0%.0f', value_now);
        else
            num_1 = sprintf('%.0f', value_now);
        end

        if value_now == 0.0
            %num_1 = sprintf('00%.0f000', value_now);
            num_1 = '0000'
        end

        num_2 = sprintf('%.0f', values_2(j)); 

        % Script names based on variables and corresponding values.
        script_names{i,j} = [var_1, num_1, '_', var_2, num_2];

    end
    
end

% Read the original file
%original_file = 'RunASE_nic5.m'; 
%file_contents = fileread(original_file);


% Loop through each gammaT value and save a new file
for i = 1:length(values_1)
    for j = 1:length(values_2)
        
        % Create the folder if it does not exist. path_scripts
        folder_name = [path_scripts, script_names{i,j}]
        
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
        end

        % Identical param file name for consistency when reading.
        full_mat = [folder_name, '/params.mat']

        % Update control values.
        ctr.meltfac = values_1(i);
        ctr.seed    = values_2(j);
        %ctr.snapshot = values_2(j);



        % We update the name of the file in accordance with the value of gammaT.
        % Avoid sign "-" in the file name as it does not compile.
        % Appropirate numering to ensure order when listing in Linux.
        
        %value_2 = ctr.snapshot;
        value_2 = ctr.seed;

        % Gamma.
        %ctr.gammaT   = values_1(i);
        %value_1 = 1e5*ctr.gammaT;
        %if value_1<10
        %    num_1 = sprintf('0%.0fe5', value_1);
        %else
        %    num_1 = sprintf('%.0fe5', value_1);
        %end

        % Melt fac.
        value_now = 1e3*values_1(i);

        if value_now<1e3
            num_1 = sprintf('0%.0f', value_now);
        else
            num_1 = sprintf('%.0f', value_now);
        end

        if value_now == 0.0
            %num_1 = sprintf('00%.0f000', value_now);
            num_1 = '0000'
        end



        num_2 = sprintf('%.0f', value_2); 


        % Input and output file names.
        name_1 = 'INIT_SSA_3'; % ASEfor_m2_80, INIT_NON_3, INIT_NON_LONG_10K
        name_2 = [var_1, num_1, '_', var_2, num_2];
    
        path_in     = [parent_path, rel_in];
        path_out    = [parent_path, rel_out, exp_1, '/', name_2, '/'];
        %path_out    = [parent_path, rel_out, exp_1, '/'];

        input       = [path_in, name_1];
        output      = [path_out, name_2];

        % Save all variables to a .mat file
        save(full_mat, 'ctr', 'input', 'output');
        
        fprintf('Saved: %s\n', full_mat);
        fprintf('Input: %s\n', input);
        fprintf('Output: %s\n', output);

    end
end