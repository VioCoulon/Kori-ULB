% Define different values for ctr.gammaT
gammaT_values = [1e-5, 2e-5, 3.e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5];
%names = ['01e5', '02e5']  

script_names = cell(length(gammaT_values), 1);


path_new_file = '/home/daniel/models/Kori-ULB/exe/thwaites/sigma_oce1/';

%num_names = [1, 2, 3, 4, 5, 6, 7, 8];
% Loop to generate filenames
for i = 1:length(gammaT_values)
    value_now = 1e5*gammaT_values(i);
    if i<10
        script_names{i} = sprintf('RunASE_lemaitre4_gamma0%.0fe5.m', value_now);
    else
        script_names{i} = sprintf('RunASE_lemaitre4_gamma%.0fe5.m', value_now);
    end
end

% Read the original file
original_file = 'RunASE_lemaitre4.m'; % Change to your actual filename
file_contents = fileread(original_file);

% Define the pattern to replace (regex for 'ctr.gammaT = <number>;')
pattern = 'ctr\.gammaT\s*=\s*[\d.eE-]+;';

% Loop through each gammaT value and save a new file
for i = 1:length(gammaT_values)
    new_value = sprintf('ctr.gammaT = %.5e;', gammaT_values(i));
    
    % Replace the old gammaT value with the new one
    new_contents = regexprep(file_contents, pattern, new_value);
    
    % Define the new filename
    %new_filename = sprintf('your_file_gammaT_%.0e.m', gammaT_values(i));
    %new_filename = script_names{i};
    full_file = [path_new_file, script_names{i}]
    
    % Save the new file
    fid = fopen(full_file, 'w');
    fwrite(fid, new_contents);
    fclose(fid);
    
    fprintf('Saved: %s\n', full_file);
end