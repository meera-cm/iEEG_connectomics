%% Generate frequency band-specific directories

addpath C:\code\spm12
% addpath C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\functions;
addpath C:\code\wjn_toolbox
spm('defaults','eeg')

%% Load table
run('ECoG_Atlas_data.m')
run('create_master_table.m')
load('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\ECoG_Atlas_Master_Table_Final.mat');

%% Assign variables

value = '_maximum_peak_amplitude'; 
freqbands = strrep(info.range_values.range_names,' ', '_'); 


% Define MRI directories for smaps and fmaps original and rh

or_sMRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\all_smaps';
or_sMRIdir = 'E:\SMAPS_RECALC\HCP_MGH_32fold_groupconnectome (Horn 2017)';
rh_sMRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';
or_fMRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\all_fmaps';
rh_fMRIdir =  'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';

%% Assign frequency bands with power spectra

% looping through the frequency bands
% for a = 1:length(freqbands)
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

%%

% ifreqnames = [freqbands {'high_beta_gt_low_beta','low_beta_gt_high_beta'}];
% disp(['N = ' num2str(sum(ifreq))])


% issue - freqnames contains low gt high beta and high gt low beta and 
% this throws an error because there aren't corresponding columns in the 
% master table. Quick fix: changing freqnames to freqbands to generate 
% original frequency bands
% for a = 1:length(ifreqnames)
for a = 1:length(freqbands)
    % create a new matrix where rows are channels and cols are freq bands
    values_matrix(:,a) = data_T.([freqbands{a} value]);
%     values_matrix(:,a) = data_T.([ifreqnames{a} value]); % e.g. data_T.alpha_maximum_peak_amplitude
    values_matrix(isnan(values_matrix(:,a)),a) = -inf;
%     values_matrix(isnan(values_matrix(:,a)),a) = -inf; % set nan vals to -infinity
end
ifreq=[];
% for a = 1:length(ifreqnames) % loop through freqband strings
for a = 1:length(freqbands)
    foi = freqbands{a};
%     foi = ifreqnames{a}; % freq of interest is a for every iteration
    fcomp = setdiff(1:length(freqbands), a);
%     fcomp = setdiff(1:length(ifreqnames),a);
    ifreq(:,a) = sum(values_matrix(:,a)>= values_matrix(:,fcomp),2) ==5;
%     ifreq(:,a)  = sum(values_matrix(:,a)>= values_matrix(:,fcomp),2)==5; % this lines literally computes the max peak amp
    subjects{a} = data_T.Patient(find(ifreq(:,a))) % selects all rows for each subject for each freqband
    %my line
%     electrodes{a} = data_T.electrodes(find(ifreq(:,a))) % theres a problem here check this
end


%% Define directories that will contain new analysis

% change path to Analysis folder
cd 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Analyses';

% create connectivity dir and 2 subdirs
mkdir connectivity str_conn
mkdir connectivity func_conn

% navigate to connectivity folder
cd connectivity

mkdir str_conn or_str_conn
mkdir str_conn rh_str_conn

mkdir func_conn or_func_conn
mkdir func_conn rh_func_conn
% should still be in conectivity folder at this point

%% Iterate through each subdir to create freqband dirs
% start the map assigning process starting with or_func_conn directory

cd func_conn
cd or_func_conn

for a = 1:length(freqbands)
    fband = freqbands(a) 
    mkdir(fband)
    fmap{:,a} = data_T.fmaps_niifiles(find(ifreq(:,a))); % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(fmap{1,a})
        filename = strcat(or_fMRIdir, '\', char(fmap{1,a}(b))) 
        copyfile(filename, char(fband), 'f');
    end 
end    
%%
cd ..

cd rh_func_conn

for a = 1:length(freqbands)
    fband = freqbands(a) 
    mkdir(fband)
    fmap{:,a} = data_T.flipped_right_fmaps(find(ifreq(:,a))); % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(fmap{1,a})
        filename = strcat(rh_fMRIdir, '\', char(fmap{1,a}(b))) 
        copyfile(filename, char(fband), 'f');
    end 
end  

%%
cd ..
cd ..
cd str_conn
cd or_str_conn

for a = 1:length(freqbands)
    fband = freqbands(a) 
    mkdir(fband)
    smap{:,a} = data_T.smaps_niifiles(find(ifreq(:,a))); % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(smap{1,a})
        filename = strcat(or_sMRIdir, '\', char(smap{1,a}(b))) 
        copyfile(filename, char(fband), 'f');
    end 
end  

%%
cd ..
cd rh_str_conn

for a = 1:length(freqbands)
    fband = freqbands(a) 
    mkdir(fband)
    smap{:,a} = data_T.flipped_right_smaps(find(ifreq(:,a))); % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(smap{1,a})
        filename = strcat(rh_sMRIdir, '\', char(smap{1,a}(b))) 
        copyfile(filename, char(fband), 'f');
    end 
end 
%% Low gt high beta original str conn

for a = 1:length(ifreqnames)
    fband = ifreqnames(a) 
    mkdir(fband)
    mkdir(ifreqnames(8))
    smap{:,8} = data_T.smaps_niifiles(find(ifreq(:,8))) % so far it's flawless
%     filename = strcat(or_sMRIdir, '\', smap
    % maybe introduce another for loop here?
    for b = 1:length(smap{1,8})
        filename = strcat(or_sMRIdir, '\', char(smap{1,8}(b))) 
        copyfile(filename, 'low_beta_gt_high_beta', 'f');
    end 
end  