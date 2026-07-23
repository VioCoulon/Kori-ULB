

% Matlab script to edit Kori run files varying certain parameter.
% Daniel Moreno Parada - March 2025.
% daniel.moreno.parada@ulb.be


% Define different values for ctr.gammaT
gammaT_values = [1e-5, 2e-5, 3.e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5];

script_names = cell(length(gammaT_values), 1);


% Experiment and variables to be changed.
exp = 'seed100';
var = 'gamma';

path_scripts = '/home/daniel/models/Kori-ULB/exe/thwaites/nic5/sigma_oce100/';
path_new_file = [path_scripts, exp, '/'];



% Loop to generate filenames
for i = 1:length(gammaT_values)
    value_now = 1e5*gammaT_values(i);
    if i<10
        num = sprintf('0%.0fe5.m', value_now);
    else
        num = sprintf('%.0fe5.m', value_now);
    end

    script_names{i} = [exp, '_', var, num];
end

% Read the original file
original_file = 'RunASE_nic5.m'; 
file_contents = fileread(original_file);

% Define the pattern to replace (regex for 'ctr.gammaT = <number>;')
pattern = 'ctr\.gammaT\s*=\s*[\d.eE-]+;';

% Loop through each gammaT value and save a new file
for i = 1:length(gammaT_values)
    new_value = sprintf('ctr.gammaT = %.5e;', gammaT_values(i));
    
    % Replace the old gammaT value with the new one
    new_contents = regexprep(file_contents, pattern, new_value);

    full_file = [path_new_file, script_names{i}]
    
    % Save the new file
    fid = fopen(full_file, 'w');
    fwrite(fid, new_contents);
    fclose(fid);
    
    fprintf('Saved: %s\n', full_file);
end