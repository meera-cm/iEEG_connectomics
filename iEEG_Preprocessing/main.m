% this is the main script, executing this will create all neccesary data.
%%
clear all
close all
clc

%%
disp('creating D');
create_D
disp('preprocessing data')
preprocess_data
disp('analysing data')
analyse_data