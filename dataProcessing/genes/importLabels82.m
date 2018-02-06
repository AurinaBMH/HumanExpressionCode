%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/genes/processedData/82parcelLabels.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2018/02/01 11:42:46

%% Import the data
[~, ~, parcelLabels] = xlsread('/Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/genes/processedData/82parcelLabels.xlsx','Sheet1');
parcelLabels(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),parcelLabels)) = {''};

