%% Creates parcellation csv for every AvgR_Fz file, calculates a correlation with an NT and calculates a correlation betweeen the correlation and beta activity

addpath C:\Code\wjn_toolbox-master
addpath C:\Code\spm12
addpath(genpath('C:\Code\leaddbs'))

%% D1
% select directory containing fmaps
fmap_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Analyses\Method_1_Correlations\fmaps';
smap_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Analyses\Method_3_Correlating_Correlations';

fmaps = dir(fullfile(fmap_dir, 'cname*.nii'))
smaps = dir(fullfile(smap_dir, 'smooth*.nii'));

parcellation_image = 'BG_Th_atlas.nii';
input_images = {fmaps.name};
% input_images = {smaps.name};

parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

fmap_tables = dir(fullfile(fmap_dir, 'cname*.csv'));
smap_tables = dir(fullfile(smap_dir, 'smooth*.csv'));
% D1_table = readtable('D1_SCH23390_hc13_kaller_HCPex (Huang 2021).csv');
DAT_table = readtable('DAT_fpcit_hc174_dukart_spect_HCPex (Huang 2021).csv');
D1_table = readtable('D1_SCH23390_hc13_kaller_HCPex (Huang 2021).csv');
D2_table = readtable('D2_fallypride_hc49_jaworska_HCPex (Huang 2021).csv');


T = table();
tempTable = table(); 
for csv = 1:length(smap_tables)
    seed_table = readtable(smap_tables(csv).name);
    [r,p] = corr(DAT_table.Value, seed_table.Value, 'type','spearman','rows', 'pairwise')
    
    tempTable.ChannelName = {smap_tables(csv).name}; 
    tempTable.rVal = r; 
    tempTable.pVal = p; 
    T = [T ; tempTable];
    writetable(T, 'smap~DAT_spearman_correlation_table.csv')
%     pause
end

%%

fmap_corr_dopa = readtable('fmap~D2_spearman_correlation_table.csv');
smap_corr_dopa = readtable('smap~DAT_spearman_correlation_table.csv');



load ECoG_Atlas_Master_Table;
% 
% size(fmap_corr_dopa.pVal)
% size(data_T.beta_maximum_peak_amplitude)

[pVal_r,pVal_p] = corr(fmap_corr_dopa.pVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman','rows', 'pairwise')
[rVal_r, rVal_p] = corr(fmap_corr_dopa.rVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman', 'rows', 'pairwise')

[pVal_r,pVal_p] = corr(smap_corr_dopa.pVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman','rows', 'pairwise')
[rVal_r, rVal_p] = corr(smap_corr_dopa.rVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman', 'rows', 'pairwise')



figure
wjn_corr_plot(fmap_corr_dopa.pVal,data_T.beta_maximum_peak_amplitude)
title('AvgRFz-D2 spearman  p values ~  beta max peak amp')
xlabel('fmap~dopa')
ylabel('beta max peak amp')
myprint('AvgRFz-D2 spearman p values ~  beta max peak amp')


figure
wjn_corr_plot(smap_corr_dopa.rVal,data_T.beta_total_integral)
title('smap-DAT spearman r values ~ beta total int')
xlabel('smap~dopa')
ylabel('beta total int')
myprint('smap-DAT spearman r values ~ beta total int')



%% Only significant r values
D2_table = readtable('fmap~D2_spearman_correlation_table.csv');
D2_table_copy = readtable('fmap~D2_spearman_correlation_table - Copy.csv');
D2_table_copy.rVal(D2_table_copy.pVal < 0.05) = nan;

D1_table_copy = readtable('fmap~D1_spearman_correlation_table_sig_vals.csv');
D1_table_copy.rVal(D1_table_copy.pVal <0.05) = nan;

DAT_table_sig = readtable('fmap~DAT_spearman_correlation_table_sig_vals.csv');
DAT_table_sig.rVal(DAT_table_sig.pVal < 0.05) = nan;

fmap_corr_dopa = readtable('fmap~DAT_spearman_correlation_table_sig_vals.csv');

[pVal_r,pVal_p] = corr(fmap_corr_dopa.pVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman','rows', 'pairwise')
[rVal_r, rVal_p] = corr(fmap_corr_dopa.rVal, data_T.beta_maximum_peak_amplitude, 'type', 'spearman', 'rows', 'pairwise')

% 
% figure
% wjn_corr_plot(fmap_corr_dopa.pVal,data_T.beta_maximum_peak_amplitude)
% title('AvgRFz-DAT spearman significant p values ~  beta max peak amp')
% xlabel('fmap~dopa')
% ylabel('beta max peak amp')
% myprint('AvgRFz-D2 spearman significant p values ~  beta max peak amp')


figure
wjn_corr_plot(fmap_corr_dopa.rVal,data_T.beta_maximum_peak_amplitude)
title('AvgRFz-DAT spearman significant r values ~  beta max peak amp')
xlabel('fmap~dopa')
ylabel('beta max peak amp')
myprint('AvgRFz-DAT spearman significant r values ~  beta max peak amp')


%%
t = splitvars(table(reshape(1:20, 5, 4)))
% set values below 8 in column 2 to nan
colN = 2; 
cutoff = 8;
t.(colN)(t.(colN) < cutoff) = nan


%% D2

fmap_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Analyses\Method_1_Correlations\MRI_Dopa_corr';
fmaps = dir(fullfile(fmap_dir, 'cname*.nii'))
parcellation_image = 'HCPex (Huang 2021).nii';
input_images = {fmaps.name};

parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

fmap_tables = dir(fullfile(fmap_dir, 'cname*.csv'));
D2_table = readtable('D2_fallypride_hc49_jaworska_HCPex (Huang 2021).csv');
DAT_table = readtable('DAT_fpcit_hc174_dukart_spect_HCPex (Huang 2021).csv');
T = table();
tempTable = table(); 
for csv = 1:length(fmap_tables)
    seed_table = readtable(fmap_tables(csv).name);
    [r,p] = corr(DAT_table.Value, seed_table.Value, 'rows', 'pairwise')

    
    tempTable.ChannelName = {fmap_tables(csv).name}; 
    tempTable.rVal = r; 
    tempTable.pVal = p; 
    T = [T ; tempTable];
    writetable(T, 'fmap~DAT_correlation_table.csv')
%     pause
end



fmap_corr_dopa = readtable('fmap~DAT_correlation_table.csv');
load ECoG_Atlas_Master_Table;

size(fmap_corr_dopa.pVal)
size(data_T.beta_maximum_peak_amplitude)

[pVal_r,pVal_p] = corr(fmap_corr_dopa.pVal, data_T.beta_maximum_peak_amplitude, 'rows', 'pairwise')
[rVal_r, rVal_p] = corr(fmap_corr_dopa.rVal, data_T.beta_maximum_peak_amplitude, 'rows', 'pairwise')


figure
wjn_corr_plot(fmap_corr_dopa.pVal,data_T.beta_maximum_peak_amplitude)
title('AvgRFz-DAT p values ~ beta max peak amp')
xlabel('fmap~dopa')
ylabel('beta max peak amp')
myprint('AvgR_Fz-DAT ~ beta max peak amp')


figure
wjn_corr_plot(fmap_corr_dopa.rVal,data_T.beta_maximum_peak_amplitude)
title('AvgRFz-DAT R values ~ beta max peak amp')
xlabel('fmap~dopa')
ylabel('beta max peak amp')
myprint('AvgR_Fz-DAT ~ beta max peak amp')

