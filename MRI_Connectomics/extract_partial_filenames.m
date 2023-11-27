%% Extract files locally using partial file name

addpath 'C:\Code\wjn_toolbox';

% Define directory that will contain all the AvgR_Fz files
% lh_AvgR_Fz_maps = 'D:\ECoG_Atlas\MRI\flipped_seeds\functional_maps\flipped_seeds_maps\left_hem_niis';
master_folder = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\master_folder';

% Define main directory containing all functional maps
% all_maps = 'D:\ECoG_Atlas\MRI\flipped_seeds\functional_maps\flipped_seeds_maps';
all_maps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\master_folder\fmaps_masked';

% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(all_maps, 'flipped_*.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, master_folder);
end

%% Extract files directly in OneDrive

% New directory with AvgR_Fz files
flipped_fmaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\flipped_fmaps';

% Directory to extract from

% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(main_dir, '*AvgR_Fz.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, flipped_fmaps);
end

%%
smooth_smaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';

% Directory to extract from
main_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps';
% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(main_dir, 'smooth*'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, smooth_smaps);
end
