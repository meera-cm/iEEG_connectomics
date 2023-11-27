% Adds freqbands to existing master table

% Load required toolboxes
addpath C:\code\spm12
addpath C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\functions;
addpath C:\code\wjn_toolbox
spm('defaults','eeg')

% Load required variables by running ECoG_Atlas_data.m file
run('ECoG_Atlas_data.m');
run('create_master_table.m')

value = '_maximum_peak_amplitude'; 
freqbands = strrep(info.range_values.range_names,' ', '_'); 
% cd to seeds file
seeds_dir = 'D:\original_seeds'
seeds_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\original_seeds';

for a = 1:length(freqbands)
    % create a new matrix where rows are channels and cols are freq bands
    values_matrix(:,a) = data_T.([freqbands{a} value]); % e.g. data_T.alpha_maximum_peak_amplitude
    values_matrix(isnan(values_matrix(:,a)),a) = -inf; % set nan vals to -infinity
end
ifreq=[];
for a = 1:length(freqbands) % loop through freqband strings
    foi = freqbands{a}; % freq of interest is a for every iteration
    fcomp = setdiff(1:length(freqbands),a); 
    ifreq(:,a)  = sum(values_matrix(:,a)>= values_matrix(:,fcomp),2)==5; % this lines literally computes the max peak amp
    subjects{a} = data_T.Patient(find(ifreq(:,a))) % selects all rows for each subject for each freqband
    %my line
%     electrodes{a} = data_T.electrodes(find(ifreq(:,a))) 
end


for a = 1:length(freqbands)
    fband = freqbands(a) 
%     mkdir fband
    fmap{:,a} = data_T.Channel_name(find(ifreq(:,a))) % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(fmap{1,a})
        filename = strcat(seeds_dir, '\', char(fmap{1,a}(b))) 
        copyfile(filename, char(fband), 'f');
    end 
end   

ls = zeros(1772,1);
for a = 1:length(data_T.Channel_name);
    [val, idx] = nanmax(table2array(data_T(a,15:20)));
    ls(a) = freqbands(idx);

end

data_T.Freqband = ls


save('ECoG_Atlas_Master_Table_New.mat')