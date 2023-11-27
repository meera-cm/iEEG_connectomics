addpath C:\code\spm12
addpath D:\ECoG_Atlas\functions
addpath C:\code\wjn_toolbox
spm('defaults','eeg')
clear all, close all

%% Define Seed folder and search for seeds with coordinates starting with a minus (pos_-) on the x axis
MRIdir = 'D:\ECoG_Atlas\MRI\seeds';
[filenames,~,fullfiles] = ffind(fullfile(MRIdir,'cname_*_pos_-*.nii'));

for a = 1:length(fullfiles)
    outputfilename = ['D:\ECoG_Atlas\MRI\flipped_seeds\flipped_' filenames{a}];
    % flipud is flip up-down
    spm_imcalc(fullfiles(a),outputfilename,'flipud(i1)')
end
