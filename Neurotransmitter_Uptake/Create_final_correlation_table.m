% Correlation table all NTs all Freqbands
% Prerequisites, parcellation image, input image - smap/fmap, wjn toolbox

addpath C:\Code\wjn_toolbox-master
addpath C:\Code\spm12
addpath(genpath('C:\Code\leaddbs'))
% 
% Atlas_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\Parcellations';
% [atlas,~,at_files] = ffind(fullfile(Atlas_dir, '*.nii'));
% 
% NT_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\NTs';
% [NT,~,nt_files] = ffind(fullfile(NT_dir, '*.nii'));
% 
% Conn_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\Connectivity_maps';
% [conn_map,~,connmaps] = ffind(fullfile(Conn_dir, '*.nii'));
% 

%%
% Atlas_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\Parcellations';
Atlas_dir = 'E:\alpha_mtx\temp\Parcellations';
[atlas,~,at_files] = ffind(fullfile(Atlas_dir, '*.nii'));

% NT_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\NTs';
NT_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\PET Data\aggregated\DoSeGluGA\NT_dir';
[NT,~,nt_files] = ffind(fullfile(NT_dir, '*.csv'));

% Conn_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Combined_atlas\Connectivity_maps';
Conn_dir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\PET Data\aggregated\DoSeGluGA\connmaps';
[conn_map,~,connmaps] = ffind(fullfile(Conn_dir, '*.csv'));

% Initialise empty table
T = table();
exampleR = nan(length(nt_files),length(connmaps)*2);
%exampleP = nan(size(exampleR));
colname = {};
rowname = {};

for nt = 1:length(nt_files)
    n = 0;
    nt_table = readtable(string(nt_files(nt)));
    for con = 1:length(connmaps)
        c = readtable(string(connmaps(con)));
        [r, p] = corr(nt_table.Value, c.Value, 'rows', 'complete');
        %[r, p] = corrcoef(nt_table.Value, c.Value, 'rows', 'complete')
       if con == 1
           ridx = 1;
           pidx = 2;
       else
           ridx = con+n;
           pidx = ridx+1;
       end
        exampleR(nt, ridx) = r;
        exampleR(nt, pidx) = p;
        if nt == 1
            [~,file,~] = fileparts(connmaps{con});
            tmp = strsplit(file,'_compound');
            colname{end+1} = [tmp{1},' R val'];
            colname{end+1} = [tmp{1},' p val'];
        end
        n = n+1;
%         t.Rda = {r,p}
%         t.Pda = p
    end
    [~,row,~] = fileparts(nt_files{nt});
    tmpr = strsplit(row,'_compound');
    rowname{end+1} = tmpr{1};
%     T = [T;t];
end
new_tab = array2table(exampleR,'VariableNames',colname','RowNames',rowname');
writetable(new_tab, 'DoSeGluGA_BG_conn_matrix.csv', 'Delimiter', ',')
disp("Done")
% 
% 
% %%
% tempdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Method_2_Results\connectivity\func_conn\or_func_conn\beta\spm_output\BG_Th'
% [t,~,d] = ffind(fullfile(tempdir, '*.csv'));
% for idx = 1:length(d)
%     path = d(idx)
%     filename = t(idx)
%     prefix = 'theta'
%     newfile = fullfile(path,[prefix, '_', filename])
%     [status, msg] = copyfile(string(filename), string(newfile))
%     
% end
% 
% for n = 1:length(d)
%     oldname = [d t(n)]
%     newname = sprintf('theta','_', t{n});
%     dos(['rename "' oldname '" "' newname '"']);
% end
% 


% fmap folder
% MRIdir = 'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\data\original_data\all_smaps';

