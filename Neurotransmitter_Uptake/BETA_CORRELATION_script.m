addpath C:\code\spm12
addpath D:\ECoG_Atlas\functions
addpath C:\code\wjn_toolbox
spm('defaults','eeg')
%% Load Table
% clear all
%%
load('D:\ECoG_Atlas\MC\ECoG_Atlas_Table.mat')
load('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\ECoG_Atlas_Table.mat')
%% Identify channels with Alpha > All other freqranges
clear BetaImages BetaFiles AlphaImages AlphaFiles

value = '_maximum_peak_amplitude';

freqbands = strrep(info.range_values.range_names,' ', '_');

for a = 1:length(freqbands)
    values_matrix(:,a) = data_T.([freqbands{a} value]); 
    values_matrix(isnan(values_matrix(:,a)),a) = -inf;
end
ifreq=[];
for a = 1:length(freqbands)
    foi = freqbands{a};
    fcomp = setdiff(1:length(freqbands),a);
    ifreq(:,a)  = sum(values_matrix(:,a)>= values_matrix(:,fcomp),2)==5;
    subjects{a} = data_T.Patient(find(ifreq(:,a)))
    
%     electrodes{a} = data_T.electrodes(find(ifreq(:,a))) % theres a problem here check this
end

% unique(subjects)
% unique(electrodes)
% add low vs high and high vs low beta
ifreq(:,length(freqbands)+1)=values_matrix(:,ci('low_beta',freqbands))<values_matrix(:,ci('high_beta',freqbands));
ifreq(:,size(ifreq,2)+1)=~ifreq(:,end);
ifreqnames = [freqbands {'high_beta_gt_low_beta','low_beta_gt_high_beta'}];

disp(['N = ' num2str(sum(ifreq))])

MRIdir = 'D:\ECoG_Atlas\MRI\maps\';

for a = 1:length(ifreqnames)
    nii_files.(ifreqnames{a}) = strcat(MRIdir,data_T.niifiles(find(ifreq(:,a))));
end

%%
% check if peak present. check for every freq band if it has the highest
% peak compared to all other freqs

%%
% tag alpha and beta images
% matlabbatch={};
% % COPY AND PASTE FROM HERE
% matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\ECoG_Atlas\MC\BETA_CORRELATION'};
% matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = nii_files.beta;
% matlabbatch{1}.spm.stats.factorial_design.cov.c =  data_T.(['beta' value])(find(ifreq(:,ci('beta',ifreqnames,1))));
% matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'MAX_PEAK_BETA';
% matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
% matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
% matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
% matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
% matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
% matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
% matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
% matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
% matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
% matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
% matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
% matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% % COPY & PASTE TILL HERE
% spm_jobman('run',matlabbatch)    

%% Plot all spectra involved in the Connectivity statistic
ibeta= find(ifreq(:,ci('beta',ifreqnames,1)));
figure
hold on
plot(F,mean(data_T.MeanPower(ibeta,:)), 'color', '#EDB120','linewidth',5)
title(['Beta N = ' num2str(numel(ibeta))]);
xlabel('Frequency [Hz]')
ylabel('Mean Power [a.u.]')
% xlim([0 80])
% ylim([0 300])
legend({'Beta'})
print('BETA_CORRELATION_power.png','-dpng')

figure
hold on
plot(F,mean(data_T.ModeledPower(ibeta,:)), 'color', '#EDB120','linewidth',5)
title(['Beta N = ' num2str(numel(ibeta))]);
xlabel('Frequency [Hz]')
ylabel('Modeled Power [a.u.]')
% xlim([0 80])
% ylim([0 300])
% legend({'Theta','Alpha','Beta','Gamma'})
legend({'Beta'})
print('BETA_CORRELATION_power.png','-dpng')



