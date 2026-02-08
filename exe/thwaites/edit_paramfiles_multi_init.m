

% Matlab script to edit Kori run files varying certain parameter.
% Daniel Moreno Parada - March 2025.
% daniel.moreno.parada@ulb.be


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORWARD RUN (CONSTANT FORCING).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First SIA inversion (As)
ctr.imax     = 472;
ctr.jmax     = 400;
ctr.delta    = 2.e3;
ctr.nsteps   = 4001; % 4001
ctr.dt       = 5;
ctr.Ao       = 5.0e-17; % 5.0e-17
ctr.m        = 3;
ctr.Asin     = zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.calving  = 2;
ctr.Tinit    = 1;    % Initialization of temperature field from semi-analytical steady-state temperature solution
ctr.Tcalc    = 2;    % Calculate temperature field and thermomechanical coupling, , i.e. A = f (T)
ctr.basin    = 1;    % Run the model for a specific basin.
ctr.inverse  = 1;  % optimization of basal sliding coefficients As for the grounded ice sheet with fixed grounding line position.
ctr.runmode  = 3;  % 1: graphics; 3: no graphics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second SSA inversion (As, MeltInv) - short run
ctr.inverse      = 2;        % Optimization of basal sliding coefficients As for the grounded ice sheet and sub-shelf melt/accretion for floating ice shelves.
ctr.meltfunc     = 1;        % Beckmann and Goosse (2003) with linear dependency on the thermal forcing, following (de Boert al., 2015). 
ctr.GroundedMelt = 0;        % necessary for basins!! Vio: 0. Frank: 1.
ctr.shelf        = 1;        % Ice shelves are considered.
ctr.SSA          = 3;        % 2 Hybrid, 3 DIVA.
ctr.nsteps       = 101;
ctr.dt           = 0.01;
ctr.Tinit        = 0;        % SET 1 FOR INITIALIZATION ! Initial temperature field read from input file or when not available kept constant at values of surface temperature
ctr.Tinv         = 10;       % Time interval between updates in the optimization scheme
ctr.TinvMelt     = 0.01;
ctr.HinvMelt     = 10;

ctr.shelftune=0.5; % Lower values make ice shelf more viscous. In ASE we strongly underestimate ice shelves velocities.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second SSA inversion (As, MeltInv) - short run
%{  %}
ctr.nsteps   = 40001; % 10001, 40001
ctr.dt       = 0.1; % 0.1
ctr.TinvMelt = 5;
ctr.Tinv     = 5;
ctr.Hinv     = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat longer run to stabilze
%{
ctr.nsteps = 40001;  % Frank: 40000
%ctr.dt     = 0.1;   % 0.1
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define different values.
%values_1 = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
values_1 = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
values_2 = [0];

%values_2 = linspace(0, 1000, 101);
%values_2 = [values_2, 2500.0];

% Create empty dictionary to populate it with exp names.
script_names = cell(length(values_1), length(values_2));

% Experiment name. Folder name.
exp_out = 'INIT_DIVA_2';   % 'deter_revert_t500', INIT_SIA, INIT_DIVA_1, INIT_DIVA_2.

% Folder with restart files. If each experiments begins with a different restart file.
exp_in = 'INIT_DIVA_1'

% Variables.
var_1   = 'INIT_shelftune';      % gamma
var_2   = 'seed';     % seed, snapshot


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
%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lyra/stoch/tau_To_140/';
%parent_path    = '/globalsc/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/stoch/tau_To_140/';                             % Stoch.


% lemaitre4.
% exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lemaitre4/stoch/
exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lemaitre4/init/';
parent_path    = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/';
rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/stoch/smb_0/tau_To_140/'; 
rel_out        = 'ice_data/eta1e7/ground_melt_0/';                              % Stoch.

%exe_path_local = '/home/daniel/models/Kori-ULB/exe/thwaites/lemaitre4/';
%parent_path    = '/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/';
%rel_in         = 'ice_data/eta1e7/ground_melt_0/';
%rel_out        = 'output/thwaites/deter/';                             % Stoch.


% Local.
%parent_path = '/home/daniel/models/Kori-ULB/';
%rel_in      = 'output/Thwaites/Initialization/eta1e7/ground_melt_0/';
%rel_out     = 'output/Thwaites/Stochastic/test/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



path_scripts  = [exe_path_local, exp_out, '/'];



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

        % Seed.
        num_2 = sprintf('%.0f', values_2(j)); 

        % tforcing.
        %if values_2(j)<1e3
        %    num_2 = sprintf('0%.0f', values_2(j));
        %else
        %    num_2 = sprintf('%.0f', values_2(j));
        %end

        %if values_2(j)<1e2
        %    num_2 = sprintf('00%.0f', values_2(j));
        %end

        %if values_2(j) == 0.0
        %    num_2 = '0000'
        %end

        % Script names based on variables and corresponding values.
        script_names{i,j} = [var_1, num_1, '_', var_2, num_2];

    end
    
end


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
        %ctr.meltfac = values_1(i);
        %ctr.seed    = values_2(j);
        %ctr.tforcing = values_2(j);

        ctr.shelftune = values_1(i);
        ctr.seed      = values_2(j);

        % We update the name of the file in accordance with the value of gammaT.
        % Avoid sign "-" in the file name as it does not compile.
        % Appropirate numering to ensure order when listing in Linux.

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


        % Seed.
        num_2 = sprintf('%.0f', values_2(j)); 

        % tforcing.
        %if values_2(j)<1e3
        %    num_2 = sprintf('0%.0f', values_2(j));
        %else
        %    num_2 = sprintf('%.0f', values_2(j));
        %end

        %if values_2(j)<1e2
        %    num_2 = sprintf('00%.0f', values_2(j));
        %end

        %if values_2(j) == 0.0
        %    num_2 = '0000'
        %end


        % Input and output file names.
        %name_1 = 'INIT_SIA';               % Forward runs with ctr.SSA=2: INIT_SSA_3. % Init: INIT_SIA, INIT_DIVA_1.
        name_2 = [var_1, num_1, '_', var_2, num_2];
    
        % If all sims start from the same restart file.
        %path_in  = [parent_path, rel_in];   

        % Consecutive initializations.
        path_in = [parent_path, rel_out, exp_in, '/', name_2, '/']; 

        % Output directory.
        path_out = [parent_path, rel_out, exp_out, '/', name_2, '/'];

        % If all sims start from the same file: name_1: [path_in, name_1]
        % If each sim start from a given restart file: name_2: [path_in, name_2]
        input_full  = [path_in, name_2];
        output_full = [path_out, name_2];

        % Save all variables to a .mat file
        save(full_mat, 'ctr', 'input_full', 'output_full');
        
        fprintf('Saved: %s\n', full_mat);
        fprintf('input_full: %s\n', input_full);
        fprintf('output_full: %s\n', output_full);

    end
end