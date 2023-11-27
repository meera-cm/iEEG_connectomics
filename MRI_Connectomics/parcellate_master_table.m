%% Add toolboxes to path
addpath C:\code\wjn_toolbox
addpath C:\code\spm12
addpath(genpath('C:\code\leaddbs'))
spm('defaults','eeg')
addpath C:\Code\wjn_toolbox-master
%% Run parcellation
clear all, close all

output_filename = 'ECoG_Atlas_MasterParc'; % Define an output filename
% BG atlas
% output_filename = 'ECoG_Atlas_MasterParc_BG';

load ECoG_Atlas_Master_Table_Final.mat % Load the table
mni_coords = [data_T.Channel_mni_X data_T.Channel_mni_Y data_T.Channel_mni_Z]; % Extract MNI coordinates from table
data_table = data_T(:,9:32); % Define the variables that should be segmented
disp(data_T.Properties.VariableNames(9:32)) % I chose these
% whole brain parcellation
% parcellation_image = '.nii'; % define the parcellation image
% compound atlas parcellation
parcellation_image = 'compound_atlas_HCPex_SUIT_ABGT.nii';
%BG atlas
% parcellation_image = 'Atlas of the Basal Ganglia and Thalamus (He 2020) - 2mm MNI152.nii';

% Run the script
wjn_mnicoord_parcellation(output_filename,mni_coords,data_table,parcellation_image);

%% Alternative Parcellation


