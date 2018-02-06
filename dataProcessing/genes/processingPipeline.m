%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-14
%% Date modified: 2017-08-29 - noise level data added
%% Date modified: 2017-11-10 - manual gene symbol naming fixed
%% Date modified: 2018-01-10 - option to update probe--> gene assingment added
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs or updates probe to gene
%   assignment (depending on options chosen)
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%%
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
clear all; 
close all; 


options.ExcludeCBandBS =  true;
options.useCUSTprobes = false;
options.updateProbes = 'reannotator'; %'Biomart', 'reannotator', 'no'; 

S1_extractData(options)

options.probeSelections = {'Variance'};
options.parcellations = {'HCP'};
options.signalThreshold = 0.5; 
options.RNAseqThreshold = 0.2; 
options.correctDistance = true; 
options.onlyMultipleProbes = false; 
options.calculateDS = true;
options.distanceCorrection = 'Euclidean'; 
options.coexpressionFor = 'all';
options.Fit = {'removeMean'};
options.distanceThreshold = 2; % first run 30, then with the final threshold 2
options.normaliseWhat = 'Lcortex';
options.normMethod = 'scaledRobustSigmoid'; 
options.percentDS =  10;
options.doNormalise = true;
options.resolution = 'ROI'; 


S2_probes(options)
% for each parcellation first run with options.distanceThreshold = 40; 
options.distanceThreshold = 40;
S3_samples2parcellation(options)
% then use the appropriate threshold
options.distanceThreshold = 2; 
S3_samples2parcellation(options)
S4_normalisation(options)


