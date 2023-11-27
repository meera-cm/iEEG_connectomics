%% clear workspace
close all;
clear all;
clc;

%% load tf data
addpath C:\code\spm12
spm('defaults','eeg')
addpath('C:\Code\fieldtrip-20210302')
addpath('C:\Code\wjn_toolbox')
% ft_defaults
%Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');
% modfied path
Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat')


%% set globals
nplots = 4;
pdim = [2 2];
ylim_raw = [-2 4];
ylim_per = [0 2];

channels  = randi([1 Dt.nchannels],1,nplots);           % plot these channels
rand_elec = randi([1 Dt.info.electrodes.n],1,nplots);   % plot these electrodes

%% average power spectrum
MeanPower = nanmean(squeeze(Dt(:,:,:,1)),3);

%% periodic and aperiodic part from fooof parameters

% fooof_curves computes the periodic and the aperiodic parts from the 
% fooof params
[L, G] = fooof_curves(Dt);

%% plot a few time frequency signals
figure;

for i = 1:nplots
    subplot(pdim(1),pdim(2),i);
    h = pcolor(Dt.time,Dt.frequencies,log10(squeeze(Dt(channels(i),:,:,1))));
    xlabel('time [s]');
    ylabel('frequency [Hz]');
    set(h, 'EdgeColor', 'none');
    title(strcat('channel ',string(channels(i))));
    colorbar();
end

%saveas(gcf,'results/fooof_visualization/1.raw_signals.png');

%% plot power all spectra for 4 electrodes
figure;

for i=1:nplots
    subplot(pdim(1),pdim(2),i);
 
    plot(Dt.frequencies,log10(MeanPower(Dt.info.electrodes.ic == rand_elec(i),:)));
    beta_patch(ylim_raw); % Add patches
    
    legend(Dt.chanlabels(find(Dt.info.electrodes.ic == rand_elec(i))),'fontsize',5);
    %ylim(ylim_raw)
    xlabel('frequency (Hz)');
    ylabel('power');
    title(Dt.info.electrodes.names(rand_elec(i)));
end

%saveas(gcf,'results/fooof_visualization/2.power_spectra.png');

%% plot only periodic components for 4 electrodes
figure;
for i=1:nplots
    subplot(pdim(1),pdim(2),i);

    plot(Dt.frequencies, log10(MeanPower(Dt.info.electrodes.ic == ...
        rand_elec(i),:)) - L(Dt.info.electrodes.ic == rand_elec(i),:));
    beta_patch(ylim_per);
    
    legend(Dt.chanlabels(find(Dt.info.electrodes.ic == rand_elec(i))),'fontsize',5);
    ylim(ylim_per)
    xlabel('frequency (Hz)');
    ylabel('power');
    title(Dt.info.electrodes.names(rand_elec(i)));
end

% saveas(gcf,'results/fooof_visualization/3.periodic_power_spectra.png');
%% plot model of periodic components of 4 electrodes
figure;
for i=1:nplots
    subplot(pdim(1),pdim(2),i);
    
    plot(Dt.frequencies, G(Dt.info.electrodes.ic == rand_elec(i),:));
    beta_patch(ylim_per);
    legend(Dt.chanlabels(find(Dt.info.electrodes.ic == rand_elec(i))),'fontsize',5);
    ylim(ylim_per)
    xlabel('frequency (Hz)');
    ylabel('power');
    title(Dt.info.electrodes.names(rand_elec(i)));
end

%saveas(gcf,'results/fooof_visualization/4.fooof_model.png');

%close all