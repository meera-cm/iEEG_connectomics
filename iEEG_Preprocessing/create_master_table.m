%% Call data and functions
addpath C:\code\spm12
% addpath D:\ECoG_Atlas\functions
addpath C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\functions;
addpath C:\code\wjn_toolbox
spm('defaults','eeg')


%% Load info data
tf_data_path = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\iEEG\tf_spm_ecog.mat';

load dataset_info_tf

%% Rearrange multidimensional matrix to a 2D matrix
data = info.range_values;
n = 0;
data_matrix = [];
for b = 1:length(data.methods)
    for a = 1:length(data.range_names)     
        n = n+1;
        field_names{n} = strrep([data.range_names{a}, '_', data.methods{b}],...
        ' ', '_'); % creates string of whats going to be in current vector;   
        data_matrix(:,n) = data.values(:,a,b);    
    end 
end     
clear data

%% Define a new table
data_T = table;
data_T.Channel_name = info.ChannelName;
data_T.Patient = info.Patient;
data_T.Hemisphere = info.Hemisphere;
data_T.Channel_type = info.ChannelType;
data_T.Channel_region = info.ChannelRegion;
data_T.Channel_mni_X = info.ChannelPosition(:,1);
data_T.Channel_mni_Y = info.ChannelPosition(:,2);
data_T.Channel_mni_Z = info.ChannelPosition(:,3);

for a = 1:length(field_names)
    data_T.(field_names{a}) = data_matrix(:,a);
end

%% Append original functional maps to data frame
% master_folder_path = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\master_folder';
% folder = master_folder_path;

all_fmaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\all_fmaps';
folder = all_fmaps;

for a = 1:length(data_T.Channel_name)
    fmaps_niifiles{a,1} = ffind(fullfile(folder, ['*' data_T.Channel_name{a} '_*_AvgR_Fz.nii']),0);
end

% Append niifiles to data_T
data_T.fmaps_niifiles = fmaps_niifiles;

% append electrdodes to data_T
% data_T.electrodes = Dt.info.electrodes.ic;

%% Append original smoothed structural maps to data frame
smap_folderpath = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\all_smaps';
smap_folderpath = 'D:\HCP_MGH_32fold_groupconnectome (Horn 2017)'
smap_folder = smap_folderpath

for a = 1:length(data_T.Channel_name)
    smaps_niifiles{a,1} = ffind(fullfile(smap_folder, ['*' data_T.Channel_name{a} '_*struc_seed.nii']),0)
end

% Append niifiles to data_T
data_T.smaps_niifiles = smaps_niifiles;

% append electrdodes to data_T
% data_T.electrodes = Dt.info.electrodes.ic

%% Append flipped functional maps to data frame
% % Note: Appending flipped images creates empty values for non flipped seeds
% % And this doesn't seem to get resolved by the following attempts
% 
% for a = 1:length(data_T.Channel_name)
%     flipped_fmaps_niifiles{a,1} = ffind(fullfile(folder, ['flipped_*' data_T.Channel_name{a} '_*AvgR_Fz.nii']),0);
% end
% 
% % Append niifiles to data_T
% data_T.flipped_fmaps_niifiles = flipped_fmaps_niifiles;

%% Append all right fmaps from master table
% 
% master_folder_path = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\master_folder';
% folder = master_folder_path;
% for a = 1:length(data_T.Channel_name)
%     flipped_right_fmaps{a,1} = ffind(fullfile(folder, ['' data_T.Channel_name{a} '_*_AvgR_Fz.nii']),0);
% end
% 
% % Append niifiles to data_T
% data_T.flipped_right_fmaps = flipped_right_fmaps;
 %% Append flipped smoothed structural maps to data frame
% 
% for a = 1:length(data_T.Channel_name)
%     flipped_smaps_niifiles{a,1} = ffind(fullfile(folder, ['flipped_*' data_T.Channel_name{a} '_*struc_seed.nii']),0);
% end
% 
% % Append niifiles to data_T
% data_T.flipped_smaps_niifiles = flipped_smaps_niifiles;

%% All right hemisphere smoothed smaps

all_right_smaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';
folder = all_right_smaps;
for a = 1:length(data_T.Channel_name)
    flipped_right_smaps{a,1} = ffind(fullfile(folder, ['smooth_*' data_T.Channel_name{a} '_*_struc_seed.nii']),0);
%       flipped_right_smaps{a,1} = ffind(fullfile(folder, [data_T.Channel_name{a} ]),0);  
end

% Append niifiles to data_T
data_T.flipped_right_smaps = flipped_right_smaps;

%% All right hemisphere fmaps

all_right_fmaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\flipped_data\all_right_fmaps';
folder = all_right_fmaps;
for a = 1:length(data_T.Channel_name)
    flipped_right_fmaps{a,1} = ffind(fullfile(folder, ['*_' data_T.Channel_name{a} '_*_AvgR_Fz.nii']),0);
%       flipped_right_smaps{a,1} = ffind(fullfile(folder, [data_T.Channel_name{a} ]),0);  
end

% Append niifiles to data_T
data_T.flipped_right_fmaps = flipped_right_fmaps;

%% Add power spectra to the table

% Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat');
% dir_path = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\iEEG';
Dt=spm_eeg_load(tf_data_path);


% average power spectrum
MeanPower = nanmean(squeeze(Dt(:,:,:,1)),3);

% periodic and aperiodic part from fooof parameters

% fooof_curves computes the periodic and the aperiodic parts from the 
% fooof params
[L, G] = fooof_curves(Dt);

data_T.MeanPower = MeanPower;
data_T.CleanedPower = MeanPower-L;
data_T.ModeledPower = G;

%% save table
F = Dt.frequencies;
save('ECoG_Atlas_Master_Table_Final','data_T','F')
