% This script separates the images for each freqband


addpath C:\code\spm12
addpath C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\functions;
addpath C:\code\wjn_toolbox
spm('defaults','eeg')
%% Load Table
clear all
%%
load('E:\ECoG_Atlas\MC\ECoG_Atlas_Master_Table.mat')

%% Identify channels with Alpha > All other freqranges
% clear BetaImages BetaFiles AlphaImages AlphaFiles

value = '_maximum_peak_amplitude';

% strings of all the names of the frequency bands
freqbands = strrep(info.range_values.range_names,' ', '_');

% looping through the frequency bands
for a = 1:length(freqbands)
    % create a new matrix where rows are channels and cols are freq bands
    values_matrix(:,a) = data_T.([freqbands{a} value]); % e.g. data_T.alpha_maximum_peak_amplitude
    values_matrix(isnan(values_matrix(:,a)),a) = -inf; % set nan vals to -infinity
end
ifreq=[];
for a = 1:length(freqbands) % loop through freqband strings
    foi = freqbands{a}; % freq of interest is a for every iteration
    fcomp = setdiff(1:length(freqbands),a); 
    ifreq(:,a)  = sum(values_matrix(:,a)>= values_matrix(:,fcomp),2)==5; % this lines literally computes the max peak amp
    subjects{a} = data_T.Patient(find(ifreq(:,a))) % selects all rows for each subject for each freqband
    %my line
%     electrodes{a} = data_T.electrodes(find(ifreq(:,a))) % theres a problem here check this
end


% unique(subjects)
% unique(electrodes)
% add low vs high and high vs low beta
% ifreq(:,length(freqbands)+1)=values_matrix(:,ci('low_beta',freqbands))<values_matrix(:,ci('high_beta',freqbands));
% ifreq(:,size(ifreq,2)+1)=~ifreq(:,end);
ifreqnames = [freqbands {'high_beta_gt_low_beta','low_beta_gt_high_beta'}];

disp(['N = ' num2str(sum(ifreq))])

% MRIdir = 'D:\ECoG_Atlas\MRI\maps\';
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';

for a = 1:length(ifreqnames)
%     nii_files.(ifreqnames{a}) = strcat(MRIdir,data_T.niifiles(find(ifreq(:,a))));
    fmaps_niifiles(ifreqnames{a}) =  data_T.fmaps_niifiles(find(ifreq(:,a)));
    % strcat concatenates strings horizontally
end
% this last for loop is the bane of my existence. Replacing it with a hack

%% Create individual folders for each freqband - ORIGINAL THETA

% create new theta dir
% mkdir theta_fmaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';


% Extract all fmaps in theta band
theta_smaps = data_T.smaps_niifiles(find(ifreq(:,1)))

% Create new dir with these files

mkdir theta_smaps;
theta_dir_path = '';

% copy theta files into dir
for a=1:length(theta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(theta_smaps(a)));
    copyfile(filename, 'theta_smaps', 'f');
end
    
    


% 1. create filepath in for loop for each theta file
% 2. copy that file into the new dir 
% copyfile(dummyfile, 'dummydir') where dummyfile is the whole filepath
% copyfile(dummyfile, 'dummydir', 'f')

%% Create individual folders for each freqband - RIGHT HEMISPHERE THETA

% create new theta dir
mkdir theta_smaps

% fmap folder
% MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';

% Extract all fmaps in theta band
rh_theta_smaps = data_T.flipped_right_smaps(find(ifreq(:,1)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(rh_theta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(rh_theta_smaps(a)))
    copyfile(filename, 'theta_smaps', 'f');
end


%% Create individual folders for each freqband - ORIGINAL ALPHA



% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';


% Extract all fmaps in theta band
alpha_smaps = data_T.smaps_niifiles(find(ifreq(:,2)))

% Create new dir with these files

mkdir alpha_smaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(alpha_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(alpha_smaps(a)));
    copyfile(filename, 'alpha_smaps', 'f');
end
    
%% Create individual folders for each freqband - RIGHT HEMISPHERE ALPHA

% create new alpha dir
mkdir alpha_smaps

% MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';

% Extract all fmaps in theta band
rh_gamma_smaps = data_T.flipped_right_smaps(find(ifreq(:,2)));

% Create new dir with these files

% copy alpha files into dir
for a=1:length(rh_gamma_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(rh_gamma_smaps(a)))
    copyfile(filename, 'alpha_smaps', 'f');
end


%% Create individual folders for each freqband - ORIGINAL BETA

% create new theta dir
mkdir beta_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_fmaps';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,3)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'beta_smaps', 'f');
end
    
%% Create individual folders for each freqband - RIGHT HEMISPHERE BETA

% create new alpha dir
mkdir beta_smaps

MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';

% Extract all fmaps in theta band
rh_beta_smaps = data_T.flipped_right_smaps(find(ifreq(:,3)))

% Create new dir with these files

% copy alpha files into dir
for a=1:length(rh_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(rh_beta_smaps(a)))
    copyfile(filename, 'beta_smaps', 'f');
end

%% Create individual folders for each freqband - RIGHT HEMISPHERE ALPHA

% create new alpha dir
mkdir alpha_smaps

MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';

% Extract all fmaps in theta band
rh_gamma_smaps = data_T.flipped_right_smaps(find(ifreq(:,2)))

% Create new dir with these files

% copy alpha files into dir
for a=1:length(rh_gamma_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(rh_gamma_smaps(a)))
    copyfile(filename, 'alpha_fmaps', 'f');
end


%% Create individual folders for each freqband - ORIGINAL GAMMA

% create new theta dir
mkdir gamma_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_fmaps';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,6)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'gamma_smaps', 'f');
end
    

%% Create individual folders for each freqband - RIGHT HEMISPHERE GAMMA

% create new alpha dir
mkdir gamma_smaps

MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_fmaps';

% Extract all fmaps in theta band
rh_gamma_smaps = data_T.flipped_right_smaps(find(ifreq(:,6)))

% Create new dir with these files

% copy alpha files into dir
for a=1:length(rh_gamma_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(rh_gamma_smaps(a)))
    copyfile(filename, 'gamma_smaps', 'f');
end

%% Low Beta Original
% mkdir low_beta_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,4)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'low_beta_smaps', 'f');
end
    

%% Low Beta rh
mkdir low_beta_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,4)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'low_beta_smaps', 'f');
end

%% High Beta Original

% mkdir high_beta_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,5)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'high_beta_smaps', 'f');
end
    

%% High Beta rh
mkdir high_beta_smaps

% fmap folder
MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\flipped_data\all_right_smaps\all_right_smaps_smoothed';


% Extract all fmaps in theta band
high_beta_smaps = data_T.smaps_niifiles(find(ifreq(:,5)))

% Create new dir with these files

% mkdir theta_fmaps;
% theta_dir_path = '';

% copy theta files into dir
for a=1:length(high_beta_smaps)
%     filename = strcat(MRIdir, thet
    filename = strcat(MRIdir,'\', char(high_beta_smaps(a)));
    copyfile(filename, 'high_beta_smaps', 'f');
end
%%
% check if peak present. check for every freq band if it has the highest
% peak compared to all other freqs

%%
% tag alpha and beta images
matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\ECoG_Atlas\MC\AlphaVsBeta'};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = nii_files.beta;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = nii_files.alpha;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',matlabbatch)    

%% Plot all spectra involved in the Connectivity statistic

figure
plot(F,mean(data_T.ModeledPower(i,:)),'linewidth',2)
hold on
plot(F,mean(data_T.ModeledPower(~i,:)),'linewidth',2)
title({['Beta>Alpha N = ' num2str(sum(i))],['Alpha>Beta N = ' num2str(sum(~i))]})
xlabel('Frequency [Hz]')
ylabel('Modeled Power [a.u.]')
legend('Beta>Alpha','Alpha>Beta')
print('Power spectrum.png','-dpng')



