% Load the dataset into MATLAB
load('dataset.mat');
load Multivariate_NTs.mat
% Specify the response variable and predictor variables
colnames = {'beta_fmap','A4B2', 'CB1', 'DA', 'GABAa','H3', 'M1', 'mGluR5', 'MU', 'NAT', 'NMDAR','5HT'};
[~, colidx] = ismember(colnames, data.Properties.VariableNames);

response = dataset(:,1); % Response variable is in the first column
response = beta
predictors = nts
predictors = dataset(:,2:end); % Predictor variables are in the remaining columns

% Fit a multivariate linear regression model
mdl = fitlm(predictors, response);

% Display the results of the regression model
disp(mdl);

% Plot the residuals to check for any patterns or trends
plotResiduals(mdl);

% Calculate the coefficients and statistics for the model
coeffs = table2array(mdl.Coefficients);
stats = table2array(mdl.stats);

% Display the coefficients and statistics
disp(coeffs);
disp(stats);

%%
T = readtable('Multivariate_Aggregate_NTs.csv')

beta = T{:,6};
nts = T{:,12:end};
save('Multivariate_NTs.mat', 'beta', 'nts')
load Multivariate_NTs.mat


data = xlsread('Multivariate_Aggregate_NTs.xlsx');
save('data.mat', 'data')

[data, header, ~] = xlsread('Multivariate_Aggregate_NTs.xlsx');

headerRow = header(1,:)
% save('data.mat', 'data', 'row_names', 'col_names');
newData = vertcat(headerRow, dataCell);  
newData = cell2mat(newData)
save('newData.mat', 'newData')
load 'newData.mat'

myCellArray = {'1', 2, '3'; 4, '5', 6};
myMatrix = cellfun(@str2double, myCellArray);

%%
% Load the dataset
% Create a sample double array
data = [1, 2, 3; 4, 5, 6; 7, 8, 9];

% Convert the array to a table
t = array2table(predictors);

% Rename the variable names
colnames = ['A4B2', 'CB1', 'DA', 'GABAa','H3', 'M1', 'mGluR5','MU','NAT','NMDAR','5HT','VAChT'];
colnames = str2double(colnames)

% Define the dependent variable and independent variables
y = response;
X = [predictors];

% Fit the linear regression model
mdl = fitlm(X, y);

% View the model summary
disp(mdl);

% Predict the dependent variable values
y_pred = predict(mdl, X);

% Plot the predicted values against the actual values
scatter(y, y_pred);
xlabel('Actual Weight');
ylabel('Predicted Weight');

    