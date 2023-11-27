clear all;
close all;
clc;
addpath wjn_toolbox-master
addpath spm
addpath lead dbs
%%
% 0. load tf data
Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');
Dt = spm_eeg_load('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\iEEG\tf_spm_ecog.mat');

%% 1. create spherical seed images (r=5mm)
template=fullfile(spm('dir'),'canonical','single_subj_T1.nii');
mni = Dt.info.ChannelPosition;
filenames = dir('data/fmri/seeds');
filenames = dir('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds\r10_smaps');

for i = 1:Dt.nchannels
    output_name = strcat('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds\r10_smaps\cname_', Dt.info.ChannelName(i), '_pos_', ...
        sprintf('%0.1f_',Dt.info.ChannelPosition(i,:)), '.nii');
    
    if isempty(filenames) || ~contains(output_name,{filenames.name})
        wjn_spherical_roi(output_name{:},mni(i,:),10,template);
    end    
end

%% Adapted analysis with increased seed radius

template=fullfile(spm('dir'),'canonical','single_subj_T1.nii');

% Replace with the seed filenames that you want to create (n = 12 for 12
% images that did not work).
outputfilenames = {'rcname_GD001Lc_11_pos_-56.0_-37.0_-2.0_.nii',...
                    'rcname_GD001Lo_3_pos_-8.0_42.0_-13.0_.nii'};
outputfilenames = {'rcname_GD005Lw_1_pos_-43.0_-2.0_-7.0_.nii','rcname_GD017Lt_2_pos_-40.0_-2.0_-4.0_.nii'};
% Copy the coordinates for each seed into the mni variable
mni = [-43 -2 -7;...
    -40 -2 -4];

for i = 1:length(outputfilenames)
        wjn_spherical_roi(outputfilenames{i},mni(i,:),10,template);
end