%% Generate csv files from PET data
addpath C:\Code\wjn_toolbox-master
addpath C:\Code\spm12
addpath(genpath('C:\Code\leaddbs'))

% Read atlas parcellation
region_info = readtable('HCPex (Huang 2021).txt');
region_image = ea_load_nii('HCPex (Huang 2021).nii');

% Parcellate input images 
% input_images = {'D2_fallypride_hc49_jaworska.nii','resliced_glanat_flair_func_seed_AvgR_Fz.nii','resliced_glanat_flair_struc_seed.nii'};
% input image for D2

fmap_channel = 'cname_GD001Lc_11_pos_-56.0_-37.0_-2.0__func_seed_AvgR_Fz.nii'
fmap_channel = 'cname_GD001Lo_4_pos_-12.0_43.0_-13.0__func_seed_AvgR_Fz.nii'

input_images = {'D2_aggregate.nii', 'DAT_aggregate.nii'};
% input_images = {fmap_channel}
parcellation_image = 'compound_atlas_HCPex_SUIT_ABGT.nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

% MasterParc = readtable('ECoG_Atlas_MasterParc_HCPex (Huang 2021).csv');
D1_table = readtable('D1_SCH23390_hc13_kaller_HCPex (Huang 2021).csv');
seed_table = readtable('cname_GD001Lc_11_pos_-56.0_-37.0_-2.0__func_seed_AvgR_Fz_HCPex (Huang 2021).csv');
[r, p] = corrcoef(D1_table.Value,seed_table.Value, 'rows', 'complete')

% write correlation r value and p value into a table
channel_wise_corr_T = table(cellstr(fmap_channel),r(1,2), p(1,2))

%%
fmap_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Analyses\correlations\fmaps'
fmaps = dir(fullfile(fmap_dir, 'cname_*'))
parcellation_image = 'HCPex (Huang 2021).nii';

for seed = length(fmaps)
    fmap_channel = fmaps(seed).name
    input_image = {fmap_channel}
    parcellation_table_files=wjn_nii_parcellate(input_image, parcellation_image);
end

fmap_tables = dir(fullfile(fmap_dir, '*.csv'))
D1_table = readtable('D1_SCH23390_hc13_kaller_HCPex (Huang 2021).csv');
for csv = length(fmap_tables)
    seed_table = fmaps(csv).name
    [r,p] = corrcoef(D1_table.Value, seed_table.Value, 'rows', 'complete')
end


%%
CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);
writetable(CorrTable,'VAchT_max_peak_amp_correlation.csv','Delimiter',',')

%% Compound atlas


