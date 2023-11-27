% Correlation script for nchannels vs DOPA

%% Generate csv files from PET data
% addpath C:\Code\wjn_toolbox-master
% addpath C:\Code\spm12
addpath(genpath('C:\Code\spm12'))
addpath(genpath('C:\Code\wjn_toolbox'))
addpath(genpath('C:\Code\leaddbs\leaddbs'))

% Read atlas parcellation
region_info = readtable('HCPex (Huang 2021).txt');
region_image = ea_load_nii('HCPex (Huang 2021).nii');

% Parcellate input images 
% input_images = {'D2_fallypride_hc49_jaworska.nii','resliced_glanat_flair_func_seed_AvgR_Fz.nii','resliced_glanat_flair_struc_seed.nii'};
% input image for D2
input_images = {'alpha_str_network_n388.nii', 'beta_str_network_n982.nii', 'beta_functional_network.nii', 'alpha_functional_network.nii'};

parcellation_image = 'Automated Anatomical Labeling 3 (Rolls 2020).nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

%% Read the new parcellation table files
D2_PET = readtable('GABAa-bz_flumazenil_hc16_norgaard_compound_atlas_HCPex_SUIT_ABGT.csv');
% fMRI = readtable('resliced_glanat_flair_func_seed_AvgR_Fz_HCPex (Huang 2021).csv');
fMRI = readtable('beta_thresholded_output_spm_compound_atlas_HCPex_SUIT_ABGT.csv');
% dMRI = readtable('resliced_glanat_flair_struc_seed_HCPex (Huang 2021).csv');
dMRI = readtable('or_alpha_structural_connectivity_HCPex (Huang 2021).csv');

%% Method 2
% load parcellated images
DAT = readtable('GABAa-bz_flumazenil_hc16_norgaard_compound_atlas_HCPex_SUIT_ABGT.csv');
smap = readtable('beta_structural_connectivity_HCPex (Huang 2021).csv');
% D2_min_D1 = readtable('D2-D1_HCPex (Huang 2021).csv');
spmT = readtable('rh_beta_spmT_0001_HCPex (Huang 2021).csv');
con1 = readtable('con_0002_HCPex (Huang 2021).csv');
thresh_output = readtable('beta_thresholded_output_spm_compound_atlas_HCPex_SUIT_ABGT.csv');

% BG parcellation
DAT = readtable('DAT_fepe2i_hc6_sasaki_compound_atlas_HCPex_SUIT_ABGT.csv');
% D2_min_D1 = readtable('D2-D1_BG_Th_atlas.csv');
spmT = readtable('rh_beta_spmT_0001_BG_Th_atlas.csv');
con1 = readtable('con_0001_compound_atlas_HCPex_SUIT_ABGT.csv');
thresh_output = readtable('rh_thresh_output_BG_Th_atlas.csv');
% Visualize parcellated data

% Genetics data
DAT_exp = readtable('Dopamine_aggregate_compound_atlas_HCPex_SUIT_ABGT.csv')
spmT = readtable('sbeta_str_network_n982_compound_atlas_HCPex_SUIT_ABGT.csv');
con1 = readtable('sor_beta_structural_connectivity_compound_atlas_HCPex_SUIT_ABGT.csv');
thresh_output = readtable('thresholded_output_spm_HCPex (Huang 2021).csv');

DA = readtable('GABAa_flumazenil_hc6_dukart_compound_atlas_HCPex_SUIT_ABGT.csv');
spmT = readtable('beta_spmT_0001_compound_atlas_HCPex_SUIT_ABGT.csv');
% spmT = readtable('beta_spmT_0001_compound_atlas_HCPex_SUIT_ABGT.csv');

figure
wjn_corr_plot(wjn_gaussianize(DA.Value((spmT.Value)>0)), wjn_gaussianize(spmT.Value((spmT.Value)>0)))
title('GABA  ~ beta spmT')
xlabel('GABA')
ylabel('beta spmT')
myprint('GABA~beta_spmT')

%% Genetics correlations
% spmT
figure
wjn_corr_plot(wjn_gaussianize(DAT.Value((con1.Value)>0)), wjn_gaussianize(con1.Value((con1.Value)>0)))
title('DAT gene spatial correlation spmT (gaussianised)')
xlabel('DAT expresssion')
ylabel('beta functional connectivity original')
myprint('gene_DAT_ORspmT_gauss')

% con1
figure
wjn_corr_plot(wjn_gaussianize(DAT_exp.Value((con1.Value)>0)), wjn_gaussianize(con1.Value((con1.Value)>0)))
title('DAT gene spatial correlation con1 (gaussianised)')
xlabel('DAT expresssion')
ylabel('beta functional connectivity original con1')
myprint('gene_DAT_ORcon1_gauss')
% thresholded
figure
wjn_corr_plot(wjn_gaussianize(DAT_exp.Value((thresh_output.Value)>0)), wjn_gaussianize(thresh_output.Value((thresh_output.Value)>0)))
title('DAT gene spatial correlation thresholded output (gaussianised)')
xlabel('DAT expresssion')
ylabel('beta functional connectivity original thresholded output')
myprint('gene_DAT_ORthresh_output_gauss')

% DAT PET vs genetics
figure
wjn_corr_plot(wjn_gaussianize(DAT_exp.Value((DAT.Value)>0)), wjn_gaussianize(DAT.Value((DAT.Value)>0)))
title('DAT gene ~ DAT PET spatial correlation')
xlabel('DAT gene expresssion')
ylabel('DAT PET expression')
myprint('DAT_gene_vs_PET')


%% smap correlations
figure
wjn_corr_plot(wjn_gaussianize(DAT.Value((thresh_output.Value)>0)), wjn_gaussianize(thresh_output.Value((thresh_output.Value)>0)))
title('sample')
xlabel('DAT Receptor binding')
ylabel('beta structural map connectivity RH (sum 1005 seeds)')
myprint('sample')

%% generate fmap correlations
figure
% wjn_corr_plot(D2_min_D1.Values,thresh_output.Value)
% wjn_corr_plot(wjn_gaussianize(DAT.Value), wjn_gaussianize(thresh_output.Value))
wjn_corr_plot(wjn_gaussianize(DAT_exp.Value((thresh_output.Value)>0)), wjn_gaussianize(thresh_output.Value((thresh_output.Value)>0)))
% wjn_corr_plot(wjn_gaussianize(thresh(table2array(thresh)>0)),wjn_gaussianize(DAT(table2array(thresh)>0)))
title('Genetics g-HCPex atlas Spatial correlation (gaussianised)')
xlabel('Genetics DAT Receptor binding')
ylabel('thresholded output connectivity')
myprint('Genetics g-HCPex atlas DAT - thresholded output')

figure
wjn_corr_plot(wjn_gaussianize(DAT.Value((con1.Value)>0)), wjn_gaussianize(con1.Value((con1.Value)>0)))
title('g-HCPex atlas Spatial correlation (gaussianised)')
xlabel('DAT Receptor binding')
ylabel('con 1 connectivity')
myprint('g-HCPex atlas DAT - con1')

figure
wjn_corr_plot(wjn_gaussianize(DAT.Value((spmT.Value)>0)), wjn_gaussianize(spmT.Value((spmT.Value)>0)))
title('g-HCPex atlas Spatial correlation (gaussianised)')
xlabel('DAT Receptor binding')
ylabel('spmT connectivity')
myprint('g-HCPex atlas DAT - spmT')


%% Nandu's for loop
atlases = {'HCPex','BG_Th'}
neurotransmitters = {'DAT_fepe2i_hc6_sasaki', 'D2-D1'}
for a=1:length(atlases)
    for n=1:length(neurotransmitters)
        nameCat = [atlases{a} '-' neuro{n}]
    end
end
%% D2 Fallypride  49 jaworska dataset
% Huang 2021 atlas
% MasterParc = readtable('ECoG_Atlas_MasterParc_HCPex (Huang 2021).csv');
% compound atlas
MasterParc = readtable('ECoG_Atlas_MasterParc_HCPex (Huang 2021).csv');
% BG atlas
MasterParc = readtable('ECoG_Atlas_MasterParc_compound_atlas_HCPex_SUIT_ABGT.csv')


D1_table = readtable('D1_SCH23390_hc13_kaller_BG_Th_atlas.csv');
% D2 compound atlas parcellation
D2_table = readtable('D2_fallypride_hc49_jaworska_BG_Th_atlas.csv');
D2_min_D1_table = readtable('D2-D1_HCPex (Huang 2021).csv');
D1_min_D2_table = readtable('D1-D2_HCPex (Huang 2021).csv');
% DAT compound parcellation
DAT_table = readtable('DAT_fepe2i_hc6_sasaki_BG_Th_atlas.csv');
FDOPA_table = readtable('FDOPA_fluorodopa_hc12_gomez_HCPex (Huang 2021).csv');
% Other Neurotransmitters
SerA_table = readtable('5HT1a_way_hc36_savli_HCPex (Huang 2021).csv');
SerB_table = readtable('5HT1b_p943_hc65_gallezot_HCPex (Huang 2021).csv');



% this is the line
[r, p] = corr(MasterParc.beta_maximum_peak_amplitude, D1_min_D2_table.Values, 'type', 'spearman', 'rows', 'pairwise')
% [r,p] = corr(d1_dat, betamax_amp, 'type', 'spearman', 'rows', 'pairwise')
% create matrix with values from all freqbands max peak amp

% MaxAmp = MasterParc(:,9:14);
% % transpose for correlation
% transMaxAmp = rows2vars(MaxAmp);
% % remove first row
% TableMaxAmp = transMaxAmp(:,2:428)

% individually for each freqband
% theta
[theta_r, theta_p] = corr(MasterParc.theta_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')
[alpha_r, alpha_p] = corr(MasterParc.alpha_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')
[beta_r, beta_p] = corr(MasterParc.beta_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')
[high_beta_r, high_beta_p] = corr(MasterParc.high_beta_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')
[low_beta_r, low_beta_p] = corr(MasterParc.low_beta_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')
[gamma_r, gamma_p] = corr(MasterParc.gamma_maximum_peak_amplitude, DAT_table.Values, 'type', 'spearman', 'rows', 'pairwise')


CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'BG_atlas_DAT~max_peak_amp_SPEARMAN_correlation.csv','Delimiter',',')


figure
wjn_corr_plot(MasterParc.beta_maximum_peak_amplitude, D1_min_D2_table.Values)
title('HCPex atlas Spatial correlation D1-D2 ~ beta max peak amp')
xlabel('beta max peak amp')
ylabel('D1-D2 receptor binding')
myprint('HCPex atlas D1-D2 - beta max peak amp')

%% D2 flb457 HC55 sandiego
D2_55_table = readtable('D2_flb457_hc55_sandiego_HCPex (Huang 2021).csv');

% this is the line
[r, p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_table.Value, 'rows', 'complete')

% create matrix with values from all freqbands max peak amp

MaxAmp = MasterParc(:,9:14);
% % transpose for correlation
% transMaxAmp = rows2vars(MaxAmp);
% % remove first row
% TableMaxAmp = transMaxAmp(:,2:428)

% individually for each freqband
% theta
[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, D2_55_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'D2_55_max_peak_amp_correlation.csv','Delimiter',',')

%% D1

D1_table = readtable('D1_SCH23390_hc13_kaller_HCPex (Huang 2021).csv');

% this is the line
% [r, p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_table.Value, 'rows', 'complete')

% create matrix with values from all freqbands max peak amp

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, D1_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'D1_max_peak_amp_correlation.csv','Delimiter',',')


%% DAT

DAT_table = readtable('DAT_fpcit_hc174_dukart_spect_HCPex (Huang 2021).csv');

% this is the line
% [r, p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_table.Value, 'rows', 'complete')

% create matrix with values from all freqbands max peak amp

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, DAT_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'DAT_max_peak_amp_correlation.csv','Delimiter',',')

figure
bar(DAT_table.Value)
figone(6,70)
set(gca,'XTick',DAT_table.Index(1:2:end),'XTickLabel',strrep(DAT_table.Name(1:2:end),'_',' '),'XTickLabelRotation',45,'FontSize',5)
title('DAT bar plot')
myprint('DAT_bar_plot')

figure
wjn_corr_plot(DAT_table.Value,MasterParc.beta_maximum_peak_amplitude)
title('Spatial correlation')
xlabel('DAT Receptor binding')
ylabel('max peak amp')
myprint('DAT beta max peak amp')

figure
bar(MasterParc.beta_maximum_peak_amplitude)
figone(6,70)
set(gca,'XTick',DAT_table.Index(1:2:end),'XTickLabel',strrep(DAT_table.Name(1:2:end),'_',' '),'XTickLabelRotation',45,'FontSize',5)
title('beta bar plot')
myprint('beta_bar_plot')

%% FDOPA

fdopa_table = readtable('FDOPA_fluorodopa_hc12_gomez_HCPex (Huang 2021).csv');

% this is the line
% [r, p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_table.Value, 'rows', 'complete')

% create matrix with values from all freqbands max peak amp

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, fdopa_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'FDOPA_max_peak_amp_correlation.csv','Delimiter',',')

%% Serotonin 5HT1a_way_hc_36_savli

HT1a_table = readtable('5HT1a_way_hc36_savli_HCPex (Huang 2021).csv');

% this is the line
% [r, p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, D2_table.Value, 'rows', 'complete')

% create matrix with values from all freqbands max peak amp

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, HT1a_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'5HT1a_max_peak_amp_correlation.csv','Delimiter',',')

%% Serotonin 5HT1b

HT1b_table = readtable('5HT1b_p943_hc65_gallezot_HCPex (Huang 2021).csv');

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, HT1b_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);

writetable(CorrTable,'5HT1b_max_peak_amp_correlation.csv','Delimiter',',')

%% Control NT Vesicular acetylcholine transporter (VAchT hc18)

VAchT_table = readtable('VAChT_feobv_hc18_aghourian_sum_HCPex (Huang 2021).csv');

MaxAmp = MasterParc(:,9:14);

[theta_r, theta_p] = corrcoef(MasterParc.theta_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')
[alpha_r, alpha_p] = corrcoef(MasterParc.alpha_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')
[beta_r, beta_p] = corrcoef(MasterParc.beta_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')
[high_beta_r, high_beta_p] = corrcoef(MasterParc.high_beta_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')
[low_beta_r, low_beta_p] = corrcoef(MasterParc.low_beta_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')
[gamma_r, gamma_p] = corrcoef(MasterParc.gamma_maximum_peak_amplitude, VAchT_table.Value, 'rows', 'complete')

CorrTable = table(theta_r, theta_p, alpha_r, alpha_p, beta_r, beta_p, high_beta_r, high_beta_p, low_beta_r, low_beta_p, gamma_r, gamma_p);
writetable(CorrTable,'VAchT_max_peak_amp_correlation.csv','Delimiter',',')
