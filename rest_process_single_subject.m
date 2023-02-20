% ********************************************************************** %
% Title: Preprocessing Script for Resting-State Data [Script 2]
% Authors: Armen Bagdasarov & Kenneth Roberts
% Institution: Duke University
% ********************************************************************** %

%% Script Description
% This script runs the preprocessing steps as described below for each
% subjects. This script should not be run, but the file paths need to be
% changed to make sure everything is saved where it should be saved. Script
% 1 (rest_loop_over_subjects.m) will call this script.

%% Start
function summary_info = rest_process_single_subject(varargin)

% Declare variables as global (variables that you can access in other functions)
global proj

% Reset the random number generator to try to make the results replicable 
% This produces somewhat consistent results for clean_rawdata functions 
% and ICA which can vary slightly by run
rng('default');

%% Import data
% Get mff file names 
% Data is in .mff format from EGI EEG system
mff_filename = fullfile(proj.data_location, ...
    proj.mff_filenames{proj.currentSub}); % Get name
EEG = pop_mffimport({mff_filename},{'code'}); % Import file
summary_info.currentId = {proj.currentId}; % Save subject ID in summary info

%% Remove outer ring of electrodes
% The outer ring is often noisy in developmental samples

% Outer layer of channels to be removed (n = 24)
outer_chans = {'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' ...
    'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' ...
    'E73' 'E81' 'E88' 'E94' 'E99' 'E107'};

% Remove outer_chans
EEG = pop_select(EEG, 'nochannel', outer_chans);

%% Downsample from 1000 to 250 Hz 
EEG = pop_resample(EEG, 250);

%% Insert eyes-open and eyes-closed markers
% The current data does not contain markers for when events begin and end
% So, need to insert them

% This following is a function from another script
[EEG, info] = create_rest_events(EEG);

% Save whether any blocks overlap (binary 0/1 = no/yes) in summary info
% This would be bad
summary_info.block_overlap = info.block_overlap;

% Save whether any block ended early (binary 0/1 = no/yes) in summary info
% This would be bad too
summary_info.block_truncate = info.block_truncate;

% Keep only the eyes-closed data (in other words remove the eyes-open data)
% After this step, data will be 240 seconds for each subject
EEG = pop_rmdat(EEG, {'rs_closed'}, [0 60] ,0);

% Save length of data as an additional check
summary_info.data_length_check_240 = EEG.xmax;

%% Filter

% Low-pass filter at 40 Hz and high-pass filter at 1 Hz
EEG = pop_eegfiltnew(EEG, 'hicutoff', 40); 
    % Because low-pass filter is applied at 40 Hz, we do not
    % have to worry too much about potential line noise at 60 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 1);

%% Remove 60 Hz line noise with CleanLine
% Because some participants still have line noise after filtering

% CleanLine
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',1,...
    'linefreqs',60,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,...
    'plotfigures',0,'scanforlines',1,'sigtype','Channels','taperbandwidth',2,...
    'tau',50,'verb',1,'winsize',4,'winstep',1);
% All are default parameters, except tau changed from 100 to 50, 
% which is recomended by the authors for continuous unepoched data

%% Reject bad channels

% Save variable with reduced (n = 105) channel locations
all_chan_locs = EEG.chanlocs;

% Remove online reference (Cz/E129), which is flat
% We will add it back in later when re-refernecing
EEG = pop_select(EEG, 'nochannel', {'E129'});

% Save variable with reduced (n = 104) channel locations
% This will be needed later when interpolating bad channels
reduced_chan_locs = EEG.chanlocs;

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off',...
    'WindowCriterion','off','BurstRejection','off','Distance','Euclidian',...
    'MaxMem', 60); 
        % All default settings
        % Only rejecting bad channels, everything else (e.g., ASR) turned off
        % MaxMem = 60gb for reproducibility (will vary based on computer RAM)

% Save which channels were bad in summary info
if isempty(EEG.chaninfo.removedchans) % If no bad chans...
    summary_info.bad_chans = {[]}; % ...then leave blank
else
    bad_chans = {EEG.chaninfo.removedchans(:).labels};
    summary_info.bad_chans = {strjoin(bad_chans)};

    % Plot bad channels to identify whether there are clusters of bad channels
    bad_chan_ind = find(ismember({reduced_chan_locs(:).labels}, bad_chans));
    figure; topoplot(bad_chan_ind, reduced_chan_locs, 'style', 'blank', ...
        'emarker', {'.','k',[],10}, 'electrodes', 'ptslabels');

    % Save bad channels plot
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
        'Inches', 'PaperSize', [10, 10])
    bad_chan_plot_path = '[INSER PATH]'; 
    bad_chan_plot_name = [proj.currentId '_go_bad_chans_plot'];
    saveas(gca, fullfile(bad_chan_plot_path, bad_chan_plot_name), 'png');
    close(gcf);
end

% Save number of bad channels in summary info
summary_info.n_bad_chans = length(EEG.chaninfo.removedchans);

%% Interpolate removed bad channels, add Cz/E129 back in, and re-reference to the average
EEG = pop_interp(EEG, reduced_chan_locs, 'spherical');
EEG.data(105,:) = zeros(1, EEG.pnts);
EEG.chanlocs = all_chan_locs;
EEG = pop_reref( EEG, []);

%% Remove large artifacts

% Artifact Subspace Reconstruction (ASR) + 
% additional removal of bad data periods

% First, save data before ASR
% ICA will be run later on the data post-ASR
% But ICA fields will be applied to the pre-ASR data
EEG_no_rej = EEG;

% ASR
% All default settings
% Most importantly the burst criterion is set conservatively to 20 
% and burst rejection is set to on (meaning remove instead of fix bad data)
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off',...
    'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7], 'MaxMem', 60); 
        % MaxMem set to 60gb for reproducibility 


% Save how many seconds of data is left after ASR in summary_info
% This will be important later for excluding participants 
% For example, if ICA was run only on 30 seconds of data because ASR cut
% out the rest, should get rid of the file (i.e., their data was probably 
% very noisy). Likely not enough data for ICA to be reliable.
summary_info.post_ASR_data_length = EEG.xmax;

% Rereference again to reset the data to zero-sum across channels
EEG = pop_reref( EEG, []);

% ********************************************************************** %

%% ICA
% Extended infomax ICA with PCA dimension reduction
% PCA dimension reduction is necessary because of the large number of 
% channels and relatively short amount of data
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', 30);

% Save ICA plot
ica_plot_path = '[INSERT PATH]'; 
ica_plot_name = [proj.currentId '_rest_ica_plot'];
pop_topoplot(EEG, 0, [1:30], 'Independent Components', 0, 'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_plot_path, ica_plot_name), 'png');
close(gcf);

%% Select independent components (ICs) related to eye or muscle only

% Automatic classification with ICLabel
EEG = pop_iclabel(EEG, 'default');

% Flag components with >= 70% of being eye or muscle
EEG = pop_icflag(EEG, [NaN NaN; 0.7 1; 0.7 1; NaN NaN; NaN NaN; ...
    NaN NaN; NaN NaN]);

% Select components with >= 70% of being eye or musscle
eye_prob = EEG.etc.ic_classification.ICLabel.classifications(:,3);
muscle_prob = EEG.etc.ic_classification.ICLabel.classifications(:,2);
eye_rej = find(eye_prob >= .70);
muscle_rej = find(muscle_prob >= .70);
eye_muscle_rej = [eye_rej; muscle_rej];
eye_muscle_rej = eye_muscle_rej';

% Save retained variance post-ICA in summary info
% This is important because if too much variance is lost, we might want to 
% consider excluding the participant from further analyses.[projected, pvar] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, EEG.icawinv, eye_muscle_rej);
summary_info.var_retained = 100-pvar;

% Plot only the removed components
ica_rej_plot_path = '[INSERT PATH]'; 
ica_rej_plot_name = [proj.currentId '_rest_removed_ics'];
figure % This line is necessary for if there is only 1 component to plot
pop_topoplot(EEG, 0, eye_muscle_rej, 'Independent Components', 0, ...
    'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_rej_plot_path, ica_rej_plot_name), 'png');
close(gcf); close(gcf); % Need to close twice if there is only 1 component to plot

%% Copy EEG ICA fields to EEG_no_rej and remove ICs with >= 70% of being eye or muscle
% Basically, back-projecting the ICA information from the ASR-reduced data 
% to the full pre-ASR data

EEG_no_rej.icawinv = EEG.icawinv;
EEG_no_rej.icasphere = EEG.icasphere;
EEG_no_rej.icaweights = EEG.icaweights;
EEG_no_rej.icachansind = EEG.icachansind;

EEG = EEG_no_rej; % Set EEG to the one with full data length, pre-ASR

% Remove components with >= 70% of being eye or muscle
EEG = pop_subcomp(EEG, eye_muscle_rej , 0);

if isempty(eye_muscle_rej) % If no ICs removed...
    summary_info.ics_removed = {[]}; % ...then leave blank
else
    % Save which components were removed in summary info
    summary_info.ics_removed = {num2str(eye_muscle_rej)};
end

% Save the number of components removed in summary info
summary_info.n_ics_removed = length(eye_muscle_rej);

%% Artifact rejection using TBT plugin

% First epoch the data into 1 second segments
EEG = eeg_regepochs(EEG, 'recurrence', 1, 'rmbase', NaN);

% 1. Abnormal values
% Simple voltage thresholding of -100/+100 uV
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -100 , 100 ,0 , 0.996, 1, 0);

% 2. Improbable data
% Based on joint probability, SD = 3 for both local and global thresholds 
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 1, 0, 0);

% Reject based on the epochs selected above
EEG = pop_TBT(EEG, EEG.reject.rejthreshE | EEG.reject.rejjpE, 10, 1, 0);
    % Do epoch interpolation on both types of artifact rejection at once
    % Criteria must be met in at least 10 channels for the epoch to be rejected

% Save how many epochs are left 
summary_info.n_epochs = EEG.trials;

%% Plot Channel Spectra

% Plot channel spectra 
spectra_time = EEG.xmax * 1000;
figure; pop_spectopo(EEG, 1, [0  spectra_time], 'EEG' , 'freqrange',...
    [2 80],'electrodes','off');

% Save channel spectra plot
spectra_plot_path = '[INSERT PATH]'; 
spectra_plot_name = [proj.currentId '_rest_chan_spectra'];
saveas(gca, fullfile(spectra_plot_path, spectra_plot_name), 'png');
close(gcf);

%% Save final preprocessed files
% EEGLAB-compatible .set format
set_path = '[INSERT PATH]'; 
set_name = [proj.currentId '_rest_complete'];
pop_saveset(EEG, 'filepath', set_path, 'filename', set_name);
