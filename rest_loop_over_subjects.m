% ********************************************************************** %
% Title: Preprocessing Script for Resting-State Data [Script 1]
% Authors: Armen Bagdasarov & Kenneth Roberts
% Institution: Duke University
% ********************************************************************** %

%% Script Description
% This script prepares MATLAB to process resting-state data. It loops through
% all subjects and calls the rest_process_single_subject.m script, which 
% contains all of the actual preprocessing steps. Only this script needs to
% be run, but files paths need to be changed in both scripts.

%% Prepare Workspace for Preprocessing

% Clear workspace and command window
clear
clc

% Start EEGLAB (startup file in MATLAB folder should have already added it to the path or add to path manually)
eeglab

% Declare variables as global (variables that you can access in other functions)
global proj

% Path of folder with raw data
% Data is in .mff format from EGI EEG system
proj.data_location = '[INSERT PATH]';

% Get mff file names
proj.mff_filenames = dir(fullfile(proj.data_location, '*.mff'));
proj.mff_filenames = { proj.mff_filenames(:).name };

%% Loop Over Subjects and Run rest_process_single_subject.m Script
% Make sure that file paths have been changed in both scripts.

for i = 1:length(proj.mff_filenames)
    proj.currentSub = i;
    proj.currentId = proj.mff_filenames{i};
    
    % Subject ID will be filename up to first space, or up to first '.'
    space_ind = strfind(proj.currentId, ' ');
    if ~isempty(space_ind)
        proj.currentId = proj.currentId(1:(space_ind(1)-1)); 
    else
        mff_ind = strfind(proj.currentId, '.mff');
        proj.currentId = proj.currentId(1:(mff_ind(1)-1));
    end

    if i == 1
        summary_info = rest_process_single_subject;
        summary_tab = struct2table(summary_info);
    else
        summary_info = rest_process_single_subject;
        summary_row = struct2table(summary_info); % One-row table
        summary_tab = vertcat(summary_tab, summary_row); % Append new row to table
    end
    
end

%% Write Summary Info to CSV Spreadsheet

proj.output_location = '[INSERT PATH]';
writetable(summary_tab, [proj.output_location filesep 'preprocessing_log.csv']);
