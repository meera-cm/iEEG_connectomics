% add dependencies
addpath(genpath('C:/Code/wjn_toolbox-master'))
addpath(genpath('C:/Code/wjn_toolbox'))
addpath(genpath('C:/Code/leaddbs'))
addpath(genpath('C:/Code/spm12'))

%% Parcellate and export conmaps
input_images = {'alpha_functional_network.nii', 'alpha_str_network_n388.nii','beta_str_network_n982.nii', 'beta_functional_network.nii'};

% parcellation_image = 'Automated Anatomical Labeling 3 (Rolls 2020).nii';
% parcellation_table_files=wjn_nii_parcellate(input_images,parcellation_image);
parcellation_image = 'rounded_combined_regions_atlas.nii';
parcellation_table_files = wjn_nii_parcellate(input_images, parcellation_image);

% beta_fNW = readtable('ECoG_Atlas_MasterParc_sum_Automated Anatomical Labeling 3 (Rolls 2020).csv');

%% Load parcellated conmaps
beta_fmap_AAL = readtable('beta_functional_network_Automated Anatomical Labeling 3 (Rolls 2020).csv');
beta_smap_AAL = readtable('beta_str_network_n982_Automated Anatomical Labeling 3 (Rolls 2020).csv');
alpha_fmap_AAL = readtable('alpha_functional_network_Automated Anatomical Labeling 3 (Rolls 2020).csv');
alpha_smap_AAL = readtable('alpha_str_network_n388_Automated Anatomical Labeling 3 (Rolls 2020).csv');
%% sort values and create tables
% functional connectivity
beta_fmap_sorted = sortrows(beta_fmap_AAL) % manually put intensity as col 1
writetable(beta_fmap_sorted, 'beta_fmap_ranked_AAL.csv');

alpha_fmap_sorted = (sortrows(alpha_fmap_AAL)) % manually put intensity as col 1
writetable(alpha_fmap_sorted, 'alpha_fmap_ranked_AAL.csv');

% structural connectivity
beta_smap_sorted = sortrows(beta_smap_AAL) % manually put intensity as col 1
writetable(beta_smap_sorted, 'beta_smap_ranked_AAL.csv');

alpha_smap_sorted = (sortrows(alpha_smap_AAL)) % manually put intensity as col 1
writetable(alpha_smap_sorted, 'alpha_smap_ranked_AAL.csv');

%% Create barplots LR average
% Beta fmap

ranked_vals = readtable('beta_fmap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);
% First, remove 'L' or 'R' from the names to obtain the 'base' names
baseNames = erase(names, ["_R", "_L"]);

% Get unique base names
uniqueNames = unique(baseNames);

% Prepare a new intensity array for averaged intensities
averageIntensity = zeros(length(uniqueNames), 1);

% Loop through each unique name
for i = 1:length(uniqueNames)
    % Find the indices of 'L' and 'R' for this 'base' name
    nameIndices = strcmp(baseNames, uniqueNames{i});
    
    % Calculate the average intensity for these indices
    averageIntensity(i) = mean(intensity(nameIndices));
end

% Now you need to sort the values in decreasing order before plotting
[sortedIntensity, sortIndices] = sort(averageIntensity, 'ascend');
sortedNames = uniqueNames(sortIndices);

% Now you can plot your sorted and averaged results
figure;
hold on;

% mybar = barh(1:length(sortedNames), sortedIntensity, 'FaceColor',[0.2706,0.5255,0.5608]);
mycolormap = viridis;
mycoloraxis = [-0.05 0.05];
% newcolors = nan(size(mybar.CData));
for i=1:numel(sortedIntensity)
    mybar(i) = bar(i, sortedIntensity(i),'FaceColor',till_value2color(sortedIntensity(i),mycolormap,mycoloraxis));
%    newcolors(i,:)=till_value2color(sortedIntensity(i),mycolormap,mycoloraxis);
%    set(mybar(i),'FaceColor',newcolors(i,:))
end
mybar.CData(2,:)=[1 0 0]
% set(mybar,'CData',newcolors,'CDataMapping','direct')

% set(mybar,'FaceColor',newcolors)
% mybar.CData = 
xticks(1:length(sortedNames));
xticklabels(strrep(sortedNames,'_',' '))
ytickangle(45)
xlabel('Region Name (AAL)');
ylabel('Average Beta Functional Intensities');
title('Average Beta Functional Intensity by Region');

%% beta fmap bar plots left vs right

ranked_vals = readtable('beta_fmap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);

suffixR = endsWith(names, 'R');
suffixL = endsWith(names, 'L');

% Separate and sort the data based on suffixes
[intensityR, sortedIndicesR] = sort(intensity(suffixR), 'ascend');
namesR = names(suffixR);
namesR = namesR(sortedIndicesR);

[intensityL, sortedIndicesL] = sort(intensity(suffixL), 'ascend');
namesL = names(suffixL);
namesL = namesL(sortedIndicesL);

% Plot for 'R' suffix
figure;
barh(1:length(intensityR), intensityR, 'FaceColor', '#80B3FF');
yticks(1:length(namesR));
yticklabels(strrep(namesR,'_',' '));
xlabel('Beta Power');
ylabel('Region Name (AAL)');
title('Beta Power by Region (R)');

% Plot for 'L' suffix
figure;
barh(1:length(intensityL), intensityL, 'FaceColor', '#77AC30');
yticks(1:length(namesL));
yticklabels(strrep(namesL,'_',' '));
xlabel('Beta Power');
ylabel('Region Name (AAL)');
title('Beta Power by Region (L)');

%% Alpha smaps LR avg

ranked_vals = readtable('alpha_fmap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);
% First, remove 'L' or 'R' from the names to obtain the 'base' names
baseNames = erase(names, ["_R", "_L"]);

% Get unique base names
uniqueNames = unique(baseNames);

% Prepare a new intensity array for averaged intensities
averageIntensity = zeros(length(uniqueNames), 1);

% Loop through each unique name
for i = 1:length(uniqueNames)
    % Find the indices of 'L' and 'R' for this 'base' name
    nameIndices = strcmp(baseNames, uniqueNames{i});
    
    % Calculate the average intensity for these indices
    averageIntensity(i) = mean(intensity(nameIndices));
end

% Now you need to sort the values in decreasing order before plotting
[sortedIntensity, sortIndices] = sort(averageIntensity, 'ascend');
sortedNames = uniqueNames(sortIndices);

figure;
hold on;

% mybar = barh(1:length(sortedNames), sortedIntensity, 'FaceColor',[0.2706,0.5255,0.5608]);
mycolormap = viridis;
mycoloraxis = [-0.05 0.05];
% newcolors = nan(size(mybar.CData));
for i=1:numel(sortedIntensity)
    mybar(i) = bar(i, sortedIntensity(i),'FaceColor',till_value2color(sortedIntensity(i),mycolormap,mycoloraxis));
%    newcolors(i,:)=till_value2color(sortedIntensity(i),mycolormap,mycoloraxis);
%    set(mybar(i),'FaceColor',newcolors(i,:))
end
mybar.CData(2,:)=[1 0 0]
% set(mybar,'CData',newcolors,'CDataMapping','direct')

% set(mybar,'FaceColor',newcolors)
% mybar.CData = 
xticks(1:length(sortedNames));
xticklabels(strrep(sortedNames,'_',' '))
ytickangle(45)
xlabel('Region Name (AAL)');
ylabel('Average alpha functional intensities');
title('Average alpha functional intensity by region');
%%
% Now you can plot your sorted and averaged results
figure;
hold on;

barh(1:length(sortedNames), sortedIntensity, 'FaceColor',[0.2706,0.5255,0.5608]);
yticks(1:length(sortedNames));
yticklabels(strrep(sortedNames,'_',' '))
% ytickangle(45)
xlabel('Average Alpha Functional Intensities');
ylabel('Region Name (AAL)');
title('Average Alpha Structural Intensity by Region');
%% Alpha fmaps L vs R
ranked_vals = readtable('alpha_fmap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);

suffixR = endsWith(names, 'R');
suffixL = endsWith(names, 'L');

% Separate and sort the data based on suffixes
[intensityR, sortedIndicesR] = sort(intensity(suffixR), 'ascend');
namesR = names(suffixR);
namesR = namesR(sortedIndicesR);

[intensityL, sortedIndicesL] = sort(intensity(suffixL), 'ascend');
namesL = names(suffixL);
namesL = namesL(sortedIndicesL);

% Plot for 'R' suffix
figure;
barh(1:length(intensityR), intensityR, 'FaceColor', '#80B3FF');
yticks(1:length(namesR));
yticklabels(strrep(namesR,'_',' '));
xlabel('Alpha functional values');
ylabel('Region Name (AAL)');
title('Alhpa fconn by Region (R)');

% Plot for 'L' suffix
figure;
barh(1:length(intensityL), intensityL, 'FaceColor', '#77AC30');
yticks(1:length(namesL));
yticklabels(strrep(namesL,'_',' '));
xlabel('Alpha funcitonal connectivity');
ylabel('Region Name (AAL)');

%% Beta smaps LR avg
ranked_vals = readtable('beta_smap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);
% First, remove 'L' or 'R' from the names to obtain the 'base' names
baseNames = erase(names, ["_R", "_L"]);

% Get unique base names
uniqueNames = unique(baseNames);

% Prepare a new intensity array for averaged intensities
averageIntensity = zeros(length(uniqueNames), 1);

% Loop through each unique name
for i = 1:length(uniqueNames)
    % Find the indices of 'L' and 'R' for this 'base' name
    nameIndices = strcmp(baseNames, uniqueNames{i});
    
    % Calculate the average intensity for these indices
    averageIntensity(i) = mean(intensity(nameIndices));
end

% Now you need to sort the values in decreasing order before plotting
[sortedIntensity, sortIndices] = sort(averageIntensity, 'ascend');
sortedNames = uniqueNames(sortIndices);

% Now you can plot your sorted and averaged results
figure;
hold on;

barh(1:length(sortedNames), sortedIntensity, 'FaceColor',[0.2706,0.5255,0.5608]);
yticks(1:length(sortedNames));
yticklabels(strrep(sortedNames,'_',' '))
% ytickangle(45)
xlabel('Average Beta Structural Intensities');
ylabel('Region Name (AAL)');
title('Average Beta Structural Intensity by Region');

%% Beta smaps L R separate

ranked_vals = readtable('beta_smap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);

suffixR = endsWith(names, 'R');
suffixL = endsWith(names, 'L');

% Separate and sort the data based on suffixes
[intensityR, sortedIndicesR] = sort(intensity(suffixR), 'ascend');
namesR = names(suffixR);
namesR = namesR(sortedIndicesR);

[intensityL, sortedIndicesL] = sort(intensity(suffixL), 'ascend');
namesL = names(suffixL);
namesL = namesL(sortedIndicesL);

% Plot for 'R' suffix
figure;
barh(1:length(intensityR), intensityR, 'FaceColor', '#80B3FF');
yticks(1:length(namesR));
yticklabels(strrep(namesR,'_',' '));
xlabel('Beta Power');
ylabel('Region Name (AAL)');
title('Beta Power by Region (R)');

% Plot for 'L' suffix
figure;
barh(1:length(intensityL), intensityL, 'FaceColor', '#77AC30');
yticks(1:length(namesL));
yticklabels(strrep(namesL,'_',' '));
xlabel('Beta structural intensities');
ylabel('Region Name (AAL)');
title('Beta structural intensities by Region (L)');
%% Alpha smaps LR avg

ranked_vals = readtable('alpha_smap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);
% First, remove 'L' or 'R' from the names to obtain the 'base' names
baseNames = erase(names, ["_R", "_L"]);

% Get unique base names
uniqueNames = unique(baseNames);

% Prepare a new intensity array for averaged intensities
averageIntensity = zeros(length(uniqueNames), 1);

% Loop through each unique name
for i = 1:length(uniqueNames)
    % Find the indices of 'L' and 'R' for this 'base' name
    nameIndices = strcmp(baseNames, uniqueNames{i});
    
    % Calculate the average intensity for these indices
    averageIntensity(i) = mean(intensity(nameIndices));
end

% Now you need to sort the values in decreasing order before plotting
[sortedIntensity, sortIndices] = sort(averageIntensity, 'ascend');
sortedNames = uniqueNames(sortIndices);

% Now you can plot your sorted and averaged results
figure;
hold on;

barh(1:length(sortedNames), sortedIntensity, 'FaceColor',[0.2706,0.5255,0.5608]);
yticks(1:length(sortedNames));
yticklabels(strrep(sortedNames,'_',' '))
% ytickangle(45)
xlabel('Average Alhpa Structural Intensities');
ylabel('Region Name (AAL)');
title('Average Alpha Structural Intensity by Region');

%% 

ranked_vals = readtable('alpha_smap_ranked_AAL.csv');
intensity = ranked_vals.Value;
names = ranked_vals.Name;

nanIndices = isnan(intensity);
intensity = intensity(~nanIndices);
names = names(~nanIndices);

suffixR = endsWith(names, 'R');
suffixL = endsWith(names, 'L');

% Separate and sort the data based on suffixes
[intensityR, sortedIndicesR] = sort(intensity(suffixR), 'ascend');
namesR = names(suffixR);
namesR = namesR(sortedIndicesR);

[intensityL, sortedIndicesL] = sort(intensity(suffixL), 'ascend');
namesL = names(suffixL);
namesL = namesL(sortedIndicesL);

% Plot for 'R' suffix
figure;
barh(1:length(intensityR), intensityR, 'FaceColor', '#80B3FF');
yticks(1:length(namesR));
yticklabels(strrep(namesR,'_',' '));
xlabel('Alpha structural intensities');
ylabel('Region Name (AAL)');
title('Alpha strucural intensity by Region (R)');

% Plot for 'L' suffix
figure;
barh(1:length(intensityL), intensityL, 'FaceColor', '#77AC30');
yticks(1:length(namesL));
yticklabels(strrep(namesL,'_',' '));
xlabel('Alpha structural intensities');
ylabel('Region Name (AAL)');
title('Alpha structural intensities by Region (L)');
%% Alpha smaps LR separate
