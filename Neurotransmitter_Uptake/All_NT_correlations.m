%% All NT correlations

% load toolboxes
addpath C:\Code\wjn_toolbox-master
addpath C:\Code\spm12

% load NT
input_images = {'A4B2_flubatine_hc30_hillmer.nii','CB1_aggregate_map.nii',...
    'Dopamine_aggregate_map.nii','GABAa_noraard_normalized.nii','H3_cban_hc8_gallezot.nii',...
    'M1_lsn_hc24_naganawa.nii','mGluR5_aggregate_map.nii','MU_aggregate_map.nii',...
    'NAT_aggregate_map.nii','NMDAR_mean.nii','Serotonin_aggregate_map.nii','VAChT_aggregate_map.nii'};

% WB Parcellation
parcellation_image = 'compound_atlas_HCPex_SUIT_ABGT.nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);
% BG Parcellation
parcellation_image = 'ABGT.nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

%% Load parcellated cvs
NT_csv = readtable();
fband_csv
% Calculate spatial correlation
