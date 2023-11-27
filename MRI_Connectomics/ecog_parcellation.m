%% Apply a new atlas to an image
% Here the atlas HCPex_SUIT_AGBT was applied to the PET scan,
% which means that the scan was masked.
% Each parcel of the new image is now the average voxel intensity
% of the corresponding parcel on the PET scan.
% 08.04.2022 - JVH

%% SPM is needed for Lead DBS, which we use to read and write niftis
addpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\spm12')
addpath(genpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\lead'))

%% load the atlas, the PET scan of interest, and define the output
nii_atlas = ea_load_nii('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Parcellations\compound_atlas_HCPex_SUIT_ABGT.nii');
nii_PET = ea_load_nii('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ALSP\PET Data\DAT_fepe2i_hc6_sasaki.nii.gz');

nii_out = nii_atlas;

%% read in the accopanying area label table and define the output
areas = readtable('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\Parcellations\compound_atlas_HCPex_SUIT_ABGT.txt')

areas_out = table({'id'},{'intensity'},{'label'});


%% set the output to zero (or let it take the values of the PET scan, when you do not want to average)
%nii_out.img = nii_PET.img(nii_atlas.img>0); %images MUST have the same dimension
nii_out.img = zeros(nii_out.dim); % set the out image to zeros

%% each parcel of the atlas, apply it to the PET scan, then average and sum it to the output
for p = 1:length(unique(nii_atlas.img))-1
    parcel = zeros(nii_out.dim);
    avg_intensity =  mean(nii_PET.img(nii_atlas.img==p));
    parcel(nii_atlas.img==p) = avg_intensity;
    nii_out.img = nii_out.img + parcel;
    
    % add to the tabulation
    areas_out = [areas_out;{p,avg_intensity,areas.Var2{p}}];
end

%% print out the output
nii_out.fname = 'compound_atlas_DAT_fepe2i_hc6_sasaki.nii';
nii_out.pinfo = [0;0;352];
ea_write_nii(nii_out);

%% print out the accompanying text file

writetable(areas_out,'compound_atlas_DAT_fepe2i_hc6_sasaki.csv','Delimiter',',')

%%
% Spatial_correlation_pipeline (Jonathan)

%% Add lead-dbs, spm and wjn_toolbox to path
addpath C:\code\wjn_toolbox
addpath C:\code\spm12
addpath(genpath('C:\code\leaddbs'))

%% Read atlas parcellation
region_info = readtable('HCPex (Huang 2021).txt');
region_image = ea_load_nii('HCPex (Huang 2021).nii');

%% Parcellate input images 
% input_images = {'D2_fallypride_hc49_jaworska.nii','resliced_glanat_flair_func_seed_AvgR_Fz.nii','resliced_glanat_flair_struc_seed.nii'};
input_images = {'D2_fallypride_hc49_jaworska.nii','beta_spmT_0001.nii','or_beta_structural_connectivity.nii'};
input_images = {'D2_fallypride_hc49_jaworska.nii','spmT_0001.nii','or_alpha_structural_connectivity.nii'};


parcellation_image = 'HCPex (Huang 2021).nii';
parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);

%% Read the new parcellation table files
D2_PET = readtable('D2_flb457_hc37_smith_BG_Th_atlas.csv');
% fMRI = readtable('resliced_glanat_flair_func_seed_AvgR_Fz_HCPex (Huang 2021).csv');
fMRI = readtable('thresholded_output_spm_HCPex (Huang 2021).csv');
% dMRI = readtable('resliced_glanat_flair_struc_seed_HCPex (Huang 2021).csv');
dMRI = readtable('or_beta_structural_connectivity_HCPex (Huang 2021).csv');


%% Visualize parcellated data

figure
bar(D2_PET.Value)
figone(6,70)
set(gca,'XTick',D2_PET.Index(1:2:end),'XTickLabel',strrep(D2_PET.Name(1:2:end),'_',' '),'XTickLabelRotation',45,'FontSize',5)
title('D2 smith receptor binding')
myprint('D2 smith receptor binding')


figure
bar(fMRI.Value)
figone(6,70)
set(gca,'XTick',D2_PET.Index(1:2:end),'XTickLabel',strrep(D2_PET.Name(1:2:end),'_',' '),'XTickLabelRotation',45,'FontSize',5)
title('fMRI connectivity')
myprint('fMRI connectivity')


figure
bar(dMRI.Value)
figone(6,70)
set(gca,'XTick',D2_PET.Index(1:2:end),'XTickLabel',strrep(D2_PET.Name(1:2:end),'_',' '),'XTickLabelRotation',45,'FontSize',5)
title('dMRI connectivity')
myprint('dMRI connectivity');


%% Calculate spatial correlation between D2 PET and fMRI

figure
wjn_corr_plot(D2_PET.Value,fMRI.Value)
title('Spatial correlation')
xlabel('D2 Receptor binding')
ylabel('fMRI connectivity')
myprint('PET D2 - fMRI')

%% Calculate spatial correlation between D2 PET and dMRI

figure
wjn_corr_plot(D2_PET.Value,dMRI.Value)
title('Spatial correlation')
xlabel('D2 Receptor binding')
ylabel('dMRI connectivity')
myprint('PET D2 - fMRI')


