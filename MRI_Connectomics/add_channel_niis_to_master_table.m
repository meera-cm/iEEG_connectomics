%% addpaths and load variables
addpath C:\Code\wjn_toolbox-master
addpath C:\Code\spm12
addpath(genpath('C:\Code\leaddbs'))
% run ECoG_Atlas_data file
% run create_master_table file 
% to load all vars

value = '_maximum_peak_amplitude';

% strings of all the names of the frequency bands
freqbands = strrep(info.range_values.range_names,' ', '_');

% looping through the frequency bands
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
%     electrodes{a} = data_T.electrodes(find(ifreq(:,a))) % theres a problem here check this
end


% Extract all fmaps in theta band
theta_seeds = data_T.channel_niis(find(ifreq(:,1)));
alpha_seeds = data_T.channel_niis(find(ifreq(:,2)));
beta_seeds = data_T.channel_niis(find(ifreq(:,3)));
gamma_seeds = data_T.channel_niis(find(ifreq(:,6)));

% Create new dir with these files
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds';

mkdir theta_seeds
mkdir alpha_seeds;
mkdir beta_seeds;
mkdir gamma_seeds

% copy theta files into dir
for a=1:length(theta_seeds)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(theta_seeds(a)));
    copyfile(filename, 'theta_seeds', 'f');
end

% copy alpha files into dir
for a=1:length(alpha_seeds)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(alpha_seeds(a)));
    copyfile(filename, 'alpha_seeds', 'f');
end

% copy beta files into dir
for a=1:length(beta_seeds)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(beta_seeds(a)));
    copyfile(filename, 'beta_seeds', 'f');
end

% copy gamma files into dir
for a=1:length(gamma_seeds)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(gamma_seeds(a)));
    copyfile(filename, 'gamma_seeds', 'f');
end

seed_niis = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds';
folder = seed_niis;

for a = 1:length(data_T.Channel_name)
    channel_niis{a,1} = ffind(fullfile(folder, ['*' data_T.Channel_name{a} '_*_.nii']),0);
end

% Append niifiles to data_T
data_T.channel_niis = channel_niis;