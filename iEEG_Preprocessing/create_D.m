% Loads the raw iEEG data and creates an SPM file with all the neccesary variables to do further analysis such as peak statistics.

%% clear workspace
close all;
clear all;
clc;

%% load data into spm format
d = dir('data/iEEG');

if ~contains('spm_ecog.mat',{d.name})
    alldata = importdata('data/iEEG/WakefulnessMatlabFile.mat');
    data = alldata.Data';
    channels = alldata.ChannelName;
    info  = rmfield(alldata,'Data');
    D=wjn_import_rawdata('spm_ecog',data,channels,info.SamplingFrequency);
    D=chantype(D,':','LFP');
    D.info = info;
    save(D)
else
    D=spm_eeg_load('data/iEEG/spm_ecog.mat');
end
%% perform time frequency analysis

if ~contains('tf_spm_ecog.mat',{d.name})
    df = 0.25;
    freqspec = 1:df:100;
    Dt=wjn_tf_multitaper(D.fullfile,freqspec);
    save(Dt)
else
    Dt=spm_eeg_load('data/iEEG/tf_spm_ecog.mat');
end

%% delete the 2second empty intervals from power
if nnz(isnan(Dt(1,1,:))) < 162 %this means the intervals havent been set to NaN yet
    for c=1:size(Dt,1)
        t = D(c,:) == 0;
        start = findstr([0 t], [0 1]) - 1;  %gives indices of beginning of groups
        stop = findstr([t 0], [1 0]);    %gives indices of end of groups
        for n=1:size(start,2)
            Dt(c,:,Dt.indsample(D.time(start(n))):Dt.indsample(D.time(stop(n))),1) = nan;
        end
    end
    save(Dt)
end

%% do fooof analysis on average power spectra
if ~isfield(Dt.info,'peaks')
    MeanPower = nanmean(squeeze(Dt(:,:,:,1)),3);

    % FOOOF settings
    settings = struct();
    f_range = [1, 100];

    % Run FOOOF across a group of power spectra
    fooof_results = fooof_group(Dt.frequencies', MeanPower', f_range, settings);

    % compute peak integrals
    for i=1:Dt.nchannels
        fooof_results(i).peak_int = fooof_results(i).gaussian_params(:,2).*fooof_results(i).gaussian_params(:,3)*sqrt(2*pi);
    end

    Dt.info.peaks = fooof_results;

    % concatenate all the structs to get one object for each peak 
    % (useful for later analysis)

    p = [];
    contact = [];

    for i = 1:Dt.nchannels
        p = [p, structfun(@transpose,Dt.info.peaks(i),'UniformOutput',false)];
        % keep track of which peaks map to which contact
        contact = [contact; i*ones(size(Dt.info.peaks(i).peak_int))];
    end

    peaks = [];
    peaks.background_params = [p.background_params]';
    peaks.peak_params = [p.peak_params]';
    peaks.gaussian_params = [p.gaussian_params]';
    peaks.peak_int = [p.peak_int]';
    peaks.error = [p.error]';
    peaks.r_squared = [p.r_squared]';
    peaks.ic = contact;

    Dt.info.peaks = peaks;
    save(Dt) 
end

%% save the unique electrodes (mapping from contacts to electrodes)
if ~isfield(Dt.info,'electrodes')
    electrodes = regexprep(Dt.chanlabels,'\d+$','');
    [unique_electrodes,ia,ic] = unique(electrodes,'stable');
    a.ic = ic;
    a.ia = ia;
    a.names = unique_electrodes;
    a.n = max(a.ic);
    Dt.info.electrodes = a;

    save(Dt)
end

%% z-scoring

if ~isfield(Dt.info.peaks,'z')
    % z-score per electrode (a)
    z.a.pp = [];
    z.a.pi = [];

    for i = 1:Dt.info.electrodes.n
        % for each electrode zscore peak params and peak int.
        z.a.pp = [z.a.pp; zscore(Dt.info.peaks.peak_params( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.electrodes.ic == i)),:))];
        z.a.pi = [z.a.pi; zscore(Dt.info.peaks.peak_int( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.electrodes.ic == i)),:))];
    end

    % z-score per patient (b)
    z.ab.pp = [];
    z.ab.pi = [];
    z.b.pp = [];
    z.b.pi = [];

    for i = 1:max(Dt.info.Patient)
        z.ab.pp = [z.ab.pp; zscore(z.a.pp( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.Patient == i)),:))];
        z.ab.pi = [z.ab.pi; zscore(z.a.pi( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.Patient == i)),:))];

        z.b.pp = [z.b.pp; zscore(Dt.info.peaks.peak_params( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.Patient == i)),:))];
        z.b.pi = [z.b.pi; zscore(Dt.info.peaks.peak_int( ...
        ismember(Dt.info.peaks.ic,find(Dt.info.Patient == i)),:))];
    end

    % z-score everything
    z.abc.pp = zscore(z.ab.pp);
    z.abc.pi = zscore(z.ab.pi);

    z.bc.pp = zscore(z.b.pp);
    z.bc.pi = zscore(z.b.pi);

    % add to info
    Dt.info.peaks.z = z';
    save(Dt)
end

%% save range values

if ~isfield(Dt.info,'range_values')
    ranges = [[4 8]; [8 12]; [13 30]; [13 20]; [20 30]; [30 100]];
    range_names = [string('theta') string('alpha') ...
        string('beta') string('low beta') string('high beta') string('gamma')];

    [values,methods] = get_stats(Dt,Dt.info.peaks.z.bc.pp,Dt.info.peaks.z.bc.pi,ranges);

    Dt.info.range_values.values = values;
    Dt.info.range_values.ranges = ranges;
    Dt.info.range_values.range_names = range_names;
    Dt.info.range_values.methods = methods;

    save(Dt)
end