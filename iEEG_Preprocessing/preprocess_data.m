% create functional and structural connectivity maps to electrode positions
%%
clear all;
close all;
clc;

%% 0. load tf data
Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');

%% 1. create spherical seed images (r=5mm)
%template=fullfile(spm('dir'),'canonical','single_subj_T1.nii');

% Extract single sunject nii from SPM template
template = fullfile('C:\Code\spm12\canonical\single_subj_T1.nii')

mni = Dt.info.ChannelPosition;
filenames = dir('data/fmri/seeds');

for i = 1:Dt.nchannels
    output_name = strcat('data/fmri/seeds/cname_', Dt.info.ChannelName(i), '_pos_', ...
        sprintf('%0.1f_',Dt.info.ChannelPosition(i,:)), '.nii');
    
    if isempty(filenames) || ~contains(output_name,{filenames.name})
        wjn_spherical_roi(output_name{:},mni(i,:),5,template);
    end    
end

%% 2. load these seeds in lead mapper and perform connectivity
% I deleted the gunzip function from the wjn toolbox because it caused an
% error in lead mapper and was over-writing the original gunzip function.
% also even for 100 subjects this takes over 2 days to run.

% need a cell with all filenames that have not been mapped yet
allfiles = wjn_subdir(fullfile(cd,'data/fmri/seeds','*.nii'));
fdonefiles = wjn_subdir(fullfile(cd,'data/fmri/maps','*cname*AvgR_Fz.nii'));
ddonefiles = wjn_subdir(fullfile(cd,'data/fmri/maps','*cname*struc_seed.nii'));
ftodofiles = {};
dtodofiles = {};

% go through all files
if length(fdonefiles) < length(allfiles)
    for i = 1:size(allfiles,1)

        if isempty(fdonefiles) || ~contains(strrep(strrep(allfiles(i),'.nii',...
                '_func_seed_AvgR.nii'),'seeds','maps'),fdonefiles)
            ftodofiles{end+1} = allfiles{i};
        end

        if isempty(ddonefiles) || ~contains(strrep(strrep(allfiles(i),'.nii',...
                '_struc_seed.nii'),'seeds','maps'),ddonefiles)
            if ~contains(allfiles{i},'not_working')
                dtodofiles{end+1} = allfiles{i};
            end
        end
    end

    if ~isempty(dtodofiles)
        lead_job_dmri(dtodofiles)
    end

    if ~isempty(ftodofiles)
        lead_job_fmri(ftodofiles)
    end
end

%% 3.a.1 smooth dmri
filenames = wjn_subdir(fullfile(cd,'data/fmri/maps','cname*struc_seed.nii'));
sfiles = wjn_subdir(fullfile(cd,'data/fmri/maps','scname*struc_seed.nii'));

for f = 1:size(filenames,1)
    if isempty(sfiles) || ~contains(strrep(filenames(f),'cname','scname'),sfiles)
        wjn_nii_smooth(filenames(f), [8 8 8]);
    end
end

%% 3.a.2 smooth fmri
filenames = wjn_subdir(fullfile(cd,'data/fmri/maps','cname*AvgR_Fz.nii'));
sfiles = wjn_subdir(fullfile(cd,'data/fmri/maps','scname*AvgR_Fz.nii'));

for f = 1:size(filenames,1)
    if isempty(sfiles) || ~contains(strrep(filenames(f),'cname','scname'),sfiles)
        wjn_nii_smooth(filenames(f), [8 8 8]);
    end
end

%% 3.b. mask fmri
filenames = wjn_subdir(fullfile(cd,'data/fmri/maps','scname*AvgR_Fz.nii'));
mfiles = wjn_subdir(fullfile(cd,'data/fmri/maps','mscname*AvgR_Fz.nii'));

for f = 1:size(filenames,1)
    if isempty(mfiles) || ~contains(strrep(filenames(f),'scname','mscname'),mfiles)
        wjn_nii_mask(filenames(f), 'data/fmri/greymatter_mask2mm.nii', 0.0001);
    end
end

%% 4. remove all unsmoothed and unmasked files
remove_files = wjn_subdir(fullfile(cd,'data/fmri/maps','cname*.nii'));
if ~isempty(remove_files)
    delete(remove_files{:})
end

remove_files = wjn_subdir(fullfile(cd,'data/fmri/maps','scname*AvgR_Fz.nii'));
if ~isempty(remove_files)
    delete(remove_files{:})
end