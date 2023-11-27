%% clear workspace
close all;
clear all;
clc;

%% globals
beta_ranges = [[13 20];[20 30]];

%% load tf data
addpath C:\code\spm12
spm('defaults','eeg')
addpath('C:\Code\fieldtrip-20210302')
addpath('C:\Code\wjn_toolbox')
% Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');
% ft_defaults
Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat')
[beta_vals, val_names] = get_stats(Dt,Dt.info.peaks.peak_params, Dt.info.peaks.peak_int,beta_ranges);
Dt.info.beta.vals = beta_vals;
Dt.info.beta.keys = val_names;

%% what are the differences in contacts on the same electrode?
% I want one value per electrode.
beta_vals_elec = zeros(2, Dt.info.electrodes.n);

for i=1:Dt.info.electrodes.n
    min_amp = nanmin(beta_vals(Dt.info.electrodes.ic == i,:,2));
    max_amp = nanmax(beta_vals(Dt.info.electrodes.ic == i,:,2));
    
    beta_vals_elec(:,i) = (max_amp-min_amp)./max_amp;
end

%% plot variances in the electrodes
figure;
bins = linspace(0,1,11);
histogram(beta_vals_elec(1,:));%,bins);
hold on
histogram(beta_vals_elec(2,:));%,bins);
xlabel('variance (%)')
ylabel('occurences')
title('beta peaks variance (n=419)')
legend('low','high')

%saveas(gcf,'results/electrode_variance/all.png');

%% plot variances in the electrodes per channel
figure;
types = unique(Dt.info.ChannelType);

for i=1:4
    subplot(2,2,i);
    idx = unique(Dt.info.electrodes.ic(strcmp(Dt.info.ChannelType,types(i))));
    histogram(beta_vals_elec(1,idx),bins,'normalization','probability');
    hold on
    histogram(beta_vals_elec(2,idx),bins,'normalization','probability');
    xlabel('variance (%)')
    ylabel('probability')
    title(sprintf('type %c (n=%i)', string(types(i)), size(idx,1)))
    legend('low','high')
end

%saveas(gcf,'results/electrode_variance/per_channel_type.png');

%close all