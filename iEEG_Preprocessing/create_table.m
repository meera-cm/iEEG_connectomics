% Create Table
clear all
close all
%% Call data and functions
addpath C:\code\spm12
addpath D:\ECoG_Atlas\functions
addpath C:\code\wjn_toolbox
spm('defaults','eeg')

%% Load info data

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

%% find the nifti files

folder = 'D:\ECoG_Atlas\MRI\maps';
for a = 1:length(data_T.Channel_name)
    niifiles{a,1} = ffind(fullfile(folder, ['mscname*' data_T.Channel_name{a} '_*AvgR_Fz.nii']),0);
end

% Append niifiles to data_T
data_T.niifiles = niifiles

% append electrdodes to data_T
data_T.electrodes = Dt.info.electrodes.ic
%% Add power spectra to the table

Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat');

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
save('ECoG_Atlas_Table','data_T','F')



