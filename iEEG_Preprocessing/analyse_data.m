% Output a correlation map between peaks in beta range vs structural and functional connectivity.

%%
clear all;
close all;
clc;

%% load tf data
addpath('C:\Code\wjn_toolbox')
%Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');
Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat')


%% load beta (and other freq) values from fooof analysis
keys = Dt.info.range_values.methods;
range_names = Dt.info.range_values.range_names;
vals = Dt.info.range_values.values;

% for now I only care about maximum peak amplitude
method = find(contains(keys,'maximum peak amplitude'));

%% search through folders with wjn_subdir and concatenate them 
% use the preprocessed files (smoothed and masked fmri)

filenames = wjn_subdir(fullfile(cd,'data/fmri/maps','mscname*AvgR_Fz.nii'));

disp('loading fmri to matlab file')
[niix, mapping, dummy] = get_niix(filenames,Dt);

%% load the betavalue and other variables in the order the niix was loaded
betavalue = vals(mapping,:,:);
all_spectra = nanmean(Dt(:,:,:,1),3); %takes non foofed vals
spectra = all_spectra(mapping,:);
subjects = Dt.info.Patient(mapping);
electrodes = Dt.info.electrodes.ic(mapping);
%% 1.a create difference images where beta > other ranges
spm eeg
disp('creating average connectivity images (fmri)')

% I need to do this twice for the whole and separated beta ranges.
%%
range_diff_comp_argmax(filenames, betavalue, range_names, ...
    'results/correlation_maps', method, spectra, Dt.frequencies, electrodes, subjects);
range_diff_comp_argmax(filenames, betavalue(:,[1:2 4:end],:), ...
    range_names([1:2 4:end]), 'results/correlation_maps', method, spectra, Dt.frequencies, electrodes, subjects);
%% 1.b create difference images of low vs high beta
range_diff_comp_betas(filenames, betavalue, range_names, 'results/correlation_maps', method, spectra, Dt.frequencies, electrodes, subjects);

%% 2. create a correlation map as before only for electrodes with beta

for fr=1:length(range_names)
    for met=method
        disp('creating correlation map')
        cormap(niix,betavalue(:,fr,met),dummy,range_names(fr),'results/correlation_maps');
    end
end

%% 3. compare results from 1 and 2


%% 4. exchanging nans with lowest value (but be aware that data will be skewed)
nnbetavalue = nan_to_min(betavalue);

for fr=1:length(range_names)
    for met=1:length(keys)
        disp('creating correlation map nn')
        cormap(niix,nnbetavalue(:,fr,met),dummy,range_names(fr),'results\correlation_maps','nn');
    end
end

%% 5. do the loom --> can beta be predicted from fmri?

parfor fr = 1:length(range_names)
    [r,p,~,~,~,~] = wjn_fox_loom(niix,betavalue(:,fr,method),'pearson');
    solutions{fr}= [r,p]
end

for fr = length(range_names)
    r = solutions{fr}(1);
    p = solutions{fr}(2);
    save(sprintf('results/correlation_maps/%s/loom.mat',range_names(fr)),'r','p');
end