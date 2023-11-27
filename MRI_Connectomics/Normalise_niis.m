% Add dependencies
addpath('C:\Code\spm12')


nii1 = 'resliced_D1_SCH23390_hc13_kaller.nii';
nii2 = 'resliced_D2_flb457_hc37_smith.nii'; 
nii3 = 'resliced_D2_raclopride_hc7_alakurtti.nii';
nii4 = 'resliced_DAT_fepe2i_hc6_sasaki.nii';
nii5 = 'resliced_DAT_fpcit_hc174_dukart_spect.nii';
nii6 = 'resliced_FDOPA_fluorodopa_hc12_gomez.nii';

vol1 = spm_vol(nii1);
vol2 = spm_vol(nii2);
vol3 = spm_vol(nii3);
vol4 = spm_vol(nii4);
vol5 = spm_vol(nii5);
vol6 = spm_vol(nii6);

data1 = spm_read_vols(vol1);
data2 = spm_read_vols(vol2);
data3 = spm_read_vols(vol3);
data4 = spm_read_vols(vol4);
data5 = spm_read_vols(vol5);
data6 = spm_read_vols(vol6);

%% Normalise values

data1_normalised = mat2gray(data1);
data2_normalised = mat2gray(data2);
data3_normalised = mat2gray(data3);
data4_normalised = mat2gray(data4);
data5_normalised = mat2gray(data5);
data6_normalised = mat2gray(data6);

% z score values
% 
% data1_normalised = zscore(data1);
% data2_normalised = zscore(data2);
% data3_normalised = zscore(data3);
% data4_normalised = zscore(data4);
% data5_normalised = zscore(data5);
% data6_normalised = zscore(data6);
% Calculate mean

mean_data = (data1_normalised + data2_normalised + data3_normalised + data4_normalised + data5_normalised + data6_normalised)/6;

output_nii = 'Dopamine_aggregate.nii';
vol_mean = vol6;
vol_mean.fname = output_nii
spm_write_vol(vol_mean, mean_data);
%% Aggregate D2 receptors
mean_d2 = (data2_normalised + data3_normalised)/2;

D2_aggregate = 'D2_Aggregate.nii';
vol_mean = vol2;
vol_mean.fname = D2_aggregate;
spm_write_vol(vol_mean, mean_d2)

%% Aggregate DAT receptors
mean_DAT = (data4_normalised + data5_normalised)/2;
DAT_aggregate = 'DAT_Aggregate.nii';
vol_mean = vol5;
vol_mean.fname = DAT_aggregate;
spm_write_vol(vol_mean, mean_DAT);

