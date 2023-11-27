% logic

% folder = FolderPath
% fstruct = dir("*AvgR_Fz.nii)")


% select folder with all functional maps for flipped seeds
% folder = 'D:\ECoG_Atlas\MRI\flipped_seeds\functional_maps\flipped_seeds_maps';

% Select all files with AvgR_Fz from directory flipped_seeds_maps
% fstruct = dir("*AvgR_Fz.nii")
% outputs a struct of 640 filnames of all the files that were flipped
addpath 'C:\Code\wjn_toolbox';

%% Moving flipped functional maps to new directory
% Define directory that will contain all the AvgR_Fz files
lh_AvgR_Fz_maps = 'D:\ECoG_Atlas\MRI\flipped_seeds\functional_maps\flipped_seeds_maps\left_hem_niis';

% Define man directory containing all functional maps
all_maps = 'D:\ECoG_Atlas\MRI\flipped_seeds\functional_maps\flipped_seeds_maps';

% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(all_maps, '*AvgR_Fz.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, lh_AvgR_Fz_maps);
end

%% TEST Moving structural maps to a new directory directly in OneDrive

% Create directory in folder and link path to that directory 'smaps'
dummy = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\fmri\flipped_seeds_maps\test_folder';

% Define man directory containing all functional maps
all_maps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\fmri\flipped_seeds_maps';

% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(all_maps, '*seed_T.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, dummy);
end
% there were only 640 fmaps so more than 200 are missing

%% smaps local

% Create directory in folder and link path to that directory 'smaps'
dummy = 'D:\ECoG_Atlas\MRI\all_original_maps\smaps';

% Define man directory containing all functional maps
all_maps = 'D:\ECoG_Atlas\MRI\all_original_maps\all_original_maps';

% Select files with AvgR_Fz suffix
[filenames,~,fullfiles] = ffind(fullfile(all_maps, '*struc_seed.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, dummy);
end
% success!

%% Separate all left (-) files on OneDrive

% Create directory in folder and link path to that directory 'smaps'
% left_smaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_original_smaps\left_smaps';
functional_maps_lh = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\functional_maps_lh';

% Define man directory containing all functional maps
% all_maps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_original_smaps';
all_fmaps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_fmaps';

% Select files with AvgR_Fz suffix. Select lh files by adding '-'
[filenames,~,fullfiles] = ffind(fullfile(all_fmaps, 'cname*_pos_-*.nii'));

% Loop through directory moving desired files to new directory
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, functional_maps_lh);
end
% success!
   
%%

%% 
mkdir flipped_left_and_right_images
movefile('flipped*.nii','flipped_left_and_right_images')

all_files = ffind('cname*.nii')
right_files = setdiff(all_files,left_files)




%% Partitioning left and right hemisphere functional maps in main data

% This is so that the maps can be calculated as done originally as well as
% flipped onto one hemisphere

% Define new directories for lh and rh
lh_maps = 'D:\ECoG_Atlas\MRI\maps\maps_lh'
rh_maps = 'D:\ECoG_Atlas\MRI\maps\maps_rh'

% Define directory containing all original maps
all_maps = 'D:\ECoG_Atlas\MRI\maps'

% Select left AvgR_Fz maps (with '-' prefix before x coordinate)
[filenames,~,fullfiles] = ffind(fullfile(all_maps, 'mscname_*_pos_-*AvgR_Fz.nii'));

% Loop through directory moving lh maps to lh folder
for a = 1:length(fullfiles)
    file_str = char(fullfiles(a));
    movefile(file_str, lh_maps);
end

% success!
%%
% Repeat steps for right hemisphere maps

% Select right AvgR_Fz maps (only rh files sans '-' remain after previous
% command)
all_maps = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_fmaps';

[filenames,~,fullfiles] = ffind(fullfile(all_maps, 'cname_*_AvgR_Fz.nii'));

% Loop through directory moving remaining AvgR_Fz maps (right hemisphere)
% to rh folder
for b = 1:length(fullfiles)
    file_str = char(fullfiles(b));
    movefile(file_str, rh_maps);
end







