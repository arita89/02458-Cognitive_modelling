clear all;
clc;
close all;

%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /home/silvia-neurorobotics/Documents/cognitive_modeling/week11/Stimulus_presentation_script/results.txt
%
% Auto-generated by MATLAB on 14-Nov-2019 10:36:43

%% Setup the Import Options
results_file_path = "/home/silvia-neurorobotics/Documents/cognitive_modeling/week11/Stimulus_presentation_script/results.txt";
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "jpgchipjpg", "Notsmiling", "VarName5"];
opts.VariableTypes = ["double", "double", "string", "categorical", "double"];
opts = setvaropts(opts, 3, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [3, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
results = readtable(results_file_path, opts);


% Clear temporary variables
clear opts


%% Clean the data (remove data below 200ms and above 2 sec)
toDelete = results.VarName5 < 200*10^(-3);
results(toDelete,:) = [];
toDelete = results.VarName5 >2;
results(toDelete,:) = [];

%% From binary to continuos

% for each subject normalize the reaction time 
% quality of the answer from not very smiling to very smiling
user_id = [1,2];
results.VarName6 = zeros(size(results,1),1);
for i = 1:1:size(user_id,2)
    dummy_user_id = results.VarName1 == user_id(i);
    m = mean(results(dummy_user_id,:).VarName5);
    results(dummy_user_id,:).VarName6 =  1-normalize(results(dummy_user_id,:).VarName5, 'range', [0 1]);
end

toChngsign = results.Notsmiling == 'Not smiling';
results(toChngsign,6).VarName6 = - results(toChngsign,:).VarName6;