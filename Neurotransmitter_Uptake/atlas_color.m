close all, clear all, clc
addpath C:\code\wjn_toolbox
addpath C:\code\spm12
addpath(genpath('C:\code\leaddbs'))
spm('defaults','eeg')
addpath C:\Code\wjn_toolbox-master


% T = readtable('ECoG Coordinates Bipolar.csv');
% T = readtable('C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics (1)\MC\scripts\ECoG_Atlas_Master_Table.mat');
load('E:\alpha_mtx\temp\ECoG_Atlas_Master_Table.mat');

% mni = [T.Var2 T.Var3 T.Var4];
mni = [data_T.Channel_mni_X data_T.Channel_mni_Y data_T.Channel_mni_Z]

ctx=export(gifti('BrainMesh_ICBM152Left_smoothed.gii'));


%% searches for closest vertex to your mni coordinate
nmni=[];
for a=1:size(mni,1)
    [mind(1,a),i(1,a)] = min(wjn_distance(ctx.vertices,[-abs(mni(a,1)) mni(a,2:3)]));
    nmni(a,:) = ctx.vertices(i(a),:);
end
%%
% t = readtable('Automated Anatomical Labeling 3 (Rolls 2020).txt','Delimiter',' ');
% roi = {'Parietal','Frontal','Postcentral_L','Precentral_L'};
% nii = ea_load_nii('Automated Anatomical Labeling 3 (Rolls 2020).nii');
% 
% i = ci(roi,t.Var2);
% ccc= wjn_erc_colormap;
% ccc = ccc([3 4 1 2 5:end],:);
% cc = repmat(ccc(5,:),length(ctx.vertices),1);
% for a = 1:length(roi)
%     i = ci(roi{a},t.Var2);
%     ix=[];
%     for b = 1:length(i)
%         ix = [ix;find(nii.img(:)==i(b))];
%     end
%     ix = unique(ix);
%     ixx = [];
%     for b = 1:length(ix)
%          [x,y,z]=ind2sub(size(nii.img),ix(b));
%         loc = wjn_cor2mni([x,y,z],nii.mat);
%         ixx = [ixx;find(wjn_distance(ctx.vertices,loc)<2)];
%     end
%        
%     cc(unique(ixx),:) = repmat(ccc(a,:),[length(unique(ixx)) 1]);
% end
%% Create grey surface locations
ccc = wjn_erc_colormap;
ccc = ccc([3 4 1 2 5:end],:);
% close all, 
figure('color','k')
% p=wjn_plot_surface(ctx,ccc(5,:));
% p=wjn_plot_surface(ctx, [0.6 0.5, 0.60]);
p=wjn_plot_surface(ctx, [0.85, 0.8, 0.85])
% p.FaceVertexCData = cc;
figone(40,40)
view(-90,15)
camlight 
hold on
% cm = colormap('jet');
% for a=1:size(nmni,1)
% beta = data_T.beta_maximum_peak_amplitude;
beta = data_T.Channel_region;
sphere_radius = 0.75;
color_code_variable = beta(~isnan(beta),1); %% Nan not allowed.
mni_input = nmni(~isnan(beta),:);
% wjn_plot_colored_spheres(mni_input,ccc(3,:),1.5)
wjn_plot_colored_spheres(mni_input,[0.4, 0.1, 0.4],1.5)
% end
material dull
% material([ka kd ks n sc]);
hold on 
camzoom(3)
set(gcf,'color','none')
myprint('channel_locations')
%%
ccc= wjn_erc_colormap;
ccc = ccc([3 4 1 2 5:end],:);
% close all, 
figure('color','w')
p=wjn_plot_surface(ctx,ccc(5,:));
% p.FaceVertexCData = cc;
figone(40,40)
view(-90,15)
camlight 
hold on
cm = colormap('jet');
% for a=1:size(nmni,1)
beta = data_T.beta_maximum_peak_amplitude;
sphere_radius = 0.99;
color_code_variable = beta(~isnan(beta),1); %% Nan not allowed.
mni_input = nmni(~isnan(beta),:);
wjn_plot_colored_spheres(mni_input,color_code_variable,sphere_radius)
% end
material dull
hold on 
camzoom(3)
set(gcf,'color','none')
% myprint('cortical_areas_spheres')
%% alpha
ccc= wjn_erc_colormap;
ccc = ccc([3 4 1 2 5:end],:);
% close all, 
figure('color','w')
p=wjn_plot_surface(ctx,ccc(5,:));
% p.FaceVertexCData = cc;
figone(40,40)
view(-90,15)
camlight 
hold on
cm = colormap('jet');
% for a=1:size(nmni,1)
alpha = data_T.alpha_maximum_peak_amplitude;
sphere_radius = 0.90;
color_code_variable = alpha(~isnan(alpha),1); %% Nan not allowed.
mni_input = nmni(~isnan(alpha),:);
wjn_plot_colored_spheres(mni_input,color_code_variable,sphere_radius)
% end
material dull

%% gamma
ccc= wjn_erc_colormap;
ccc = ccc([3 4 1 2 5:end],:);
% close all, 
figure('color','w')
p=wjn_plot_surface(ctx,ccc(5,:));
% p.FaceVertexCData = cc;
figone(40,40)
view(-90,15)
camlight 
hold on
cm = colormap('jet');
% for a=1:size(nmni,1)
gamma = data_T.gamma_maximum_peak_amplitude;
sphere_radius = 0.90;
color_code_variable = gamma(~isnan(gamma),1); %% Nan not allowed.
mni_input = nmni(~isnan(gamma),:);
wjn_plot_colored_spheres(mni_input,color_code_variable,sphere_radius)
% end
material dull
hold on 
camzoom(3)
set(gcf,'color','none')
hold on 
camzoom(3)
set(gcf,'color','none')

%% theta
ccc= wjn_erc_colormap;
ccc = ccc([3 4 1 2 5:end],:);
close all, 
figure('color','w')
p=wjn_plot_surface(ctx,ccc(5,:));
% p.FaceVertexCData = cc;
figone(40,40)
view(-90,15)
camlight 
hold on
cm = colormap('jet');
% for a=1:size(nmni,1)
theta = data_T.theta_maximum_peak_amplitude;
sphere_radius = 0.90;
color_code_variable = theta(~isnan(theta),1); %% Nan not allowed.
mni_input = nmni(~isnan(theta),:);
wjn_plot_colored_spheres(mni_input,color_code_variable,sphere_radius)
% end
material dull
hold on 
camzoom(3)
set(gcf,'color','none')
hold on 
camzoom(3)
set(gcf,'color','none')

