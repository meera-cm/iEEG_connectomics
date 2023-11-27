% Load libraries
addpath C:\code\wjn_toolbox;
addpath C:\code\spm12;
addpath(genpath('C:\code\leaddbs'));
spm('defaults','eeg');
addpath C:\Code\wjn_toolbox-master;
addpath C:\Code\wjn_toolbox

% Parcellate ECoG master table
output_filename = 'ECoG_Atlas_MasterParc';

load ECoG_Atlas_Master_Table.mat; % Load the table
mni_coords = [data_T.Channel_mni_X data_T.Channel_mni_Y data_T.Channel_mni_Z]; % Extract MNI coordinates from table
data_table = data_T(:,9:32); % Define the variables that should be segmented
disp(data_T.Properties.VariableNames(9:32)); % I chose these

parcellation_image = 'compound_atlas_HCPex_SUIT_ABGT.nii';
wjn_mnicoord_parcellation(output_filename,mni_coords,data_table,parcellation_image);

% Parcellate DAT and beta fmap and master table to compound atlas

input_images = {'theta_thresholded.nii', 'alpha_thresholded.nii', 'gamma_thresholded.nii'};
input_images = {'beta_thresholded_output_spm.nii'};

parcellation_image = 'compound_atlas_HCPex_SUIT_ABGT.nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

%% Structural connectivity HCPex

FDOPA = readtable('FDOPA_fluorodopa_hc12_gomez_compound_atlas_HCPex_SUIT_ABGT.csv');
smap = readtable('sor_alpha_structural_connectivity_compound_atlas_HCPex_SUIT_ABGT.csv');
figure
wjn_corr_plot(wjn_gaussianize(FDOPA.Value((smap.Value)>0)), wjn_gaussianize(smap.Value((smap.Value)>0)),'structure')
title('Beta DAT dukart ~ structural connectivity spatial correlation')
xlabel('DAT dukart Receptor binding')
ylabel('Beta structural connectivity map')
myprint('DAT_dukart_HCPex_str_conn')


% Plot by group