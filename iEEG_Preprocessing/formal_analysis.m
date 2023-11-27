% Formal description of the EEG analysis. Used to get a first intuition on
% how useful the data is before adding the MRI connectivity component.

%% clear workspace
close all;
clear all;
clc;

%% load tf data
% Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');

Dt=spm_eeg_load('D:\ECoG_Atlas\iEEG\iEEG\tf_spm_ecog.mat');

% specify which z-scores you want to use
zmode = 'abc';

pp_orig = Dt.info.peaks.peak_params;
pi_orig = Dt.info.peaks.peak_int;

if strcmp(zmode,'abc')
    pp_z = Dt.info.peaks.z.abc.pp;
    pi_z = Dt.info.peaks.z.abc.pi;
elseif strcmp(zmode,'bc')
    pp_z = Dt.info.peaks.z.bc.pp;
    pi_z = Dt.info.peaks.z.bc.pi;
end

%% 1. peak histogram
% Create a histogram with all peak frequencies regardless of the frequency
% band. Calculate a linear correlation between number of peaks and peak
% frequency to see whether there is a trend towards more peaks in low or
% high frequencies.

figure('units','normalized','outerposition',[0 0 1 1]);

pf = Dt.info.peaks.peak_params(:,1);

beta_patch([0 250]);
h = histogram(pf,max(Dt.frequencies),'FaceColor','black');%,'FaceAlpha',0.2);
frequencies = h.BinEdges(1:end-1) + h.BinWidth/2;
R = corrcoef(frequencies,h.Values);

title(sprintf('Number of Peaks Decrease with Peak Frequency (R = %.3f)', R(1,2)));
xlabel('peak frequency');
ylabel('ocurrences');
legend('low beta','high beta');

% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/1.png'));

%% plot average peak amplitude for each frequency
figure;

freqs = 1:max(Dt.frequencies);
amps = zeros(size(freqs,2),1);

for f = freqs
    amp = Dt.info.peaks.peak_params((Dt.info.peaks.peak_params(:,1)>=f-0.5 & ...
        Dt.info.peaks.peak_params(:,1)<f+0.5),2);
    if isempty(amp)
        amps(f) = NaN;
    else
        amps(f) = mean(amp);
    end
end

xlabel('frequency');
ylabel('mean peak amplitude');
beta_patch([0 0.9]);
plot(freqs,amps);

% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/1b.png'));

%% 2. Create normalized values 
% and see potential effects

range_names = [string('low beta') string('high beta')];
ranges = [[13 20];[20 30]];

% plot stats for regular values
[vals,keys] = get_stats(Dt,pp_orig,pi_orig,ranges);
plot_stat_hists(vals,keys,range_names);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/2_hist_orig.png'));
plot_stat_cor(vals,keys,range_names);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/2_cor_orig.png'));

% plot stats for z-scores
[vals,keys] = get_stats(Dt,pp_z,pi_z,ranges);
plot_stat_hists(vals,keys,range_names);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/2_hist_z.png'));
plot_stat_cor(vals,keys,range_names);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/2_cor_z.png'));

%% 3. Recreate the spatial distribution plots with these normalized values
% for theta 4-8, alpha 8-12 , low beta 13-20 and high beta 20-30 hz bands.

range_names = [string('theta') string('alpha') string('low beta') string('high beta')];
ranges = [[4 8];[8 12];[13 20];[20 30]];
[vals,keys] = get_stats(Dt,pp_z,pi_z,ranges);

plot_beta_on_brain(Dt,vals,keys,ranges,range_names,[ ... %string('number of peaks') ...
    string('maximum peak amplitude') string('integral of highest peak') ...
    string('total integral')],true);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/3.png'));

%% 4. Recreate the spatial distribution plots for median split values.
% % 
val_medsplit = vals;
val_medsplit(vals < nanmedian(vals)) = 0;
val_medsplit(vals >= nanmedian(vals)) = 1;

% replace NaNs with zeros as no peaks are found there, so plotted are the
% regions where peaks are found with amplitudes higher than the median.
val_medsplit(isnan(vals)) = 0;

plot_beta_on_brain(Dt,val_medsplit,keys,ranges,range_names,[ ... %string('number of peaks') ...
    string('maximum peak amplitude') string('integral of highest peak') ...
    string('total integral')],false);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/4_medsplit.png'));

%% 5. Correlate the values for each band with x, y and z coordinates and
% distance to center to see general trends inspatial distribution.

plot_cor_xyz(vals,keys,range_names,Dt.info.ChannelPosition);
saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/5_xyz.png'));
plot_cor_cdist(vals,keys,range_names,Dt.info.ChannelPosition);
% saveas(gcf,strcat(strcat('results/formal_analysis/',zmode),'/5_cdist.png'));

%%
close all