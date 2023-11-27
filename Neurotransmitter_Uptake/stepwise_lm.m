
%% Create fband + dopa table



table = readtable('stepwise_lm_fconn.xlsx')
% model1 = stepwiselm(table, 'constant', 'ResponseVar', 'DA_aggregate')

model1 = stepwiselm(table, 'constant', 'ResponseVar', 'DA_aggregate',...
    'PredictorVars', {'theta_fconn', 'alpha_fconn', 'beta_fconn',...
    'beta_sconn', 'alpha_sconn', 'beta_sconn'});

%% GABA

table = readtable('stepwise_lm_GABA.xlsx');
model2 = stepwiselm(table, 'constant', 'ResponseVar', 'GABA', 'PredictorVars',...
    {'theta_fconn', 'alpha_fconn', 'beta_fconn',...
    'beta_sconn', 'alpha_sconn', 'beta_sconn'});
    

model2
