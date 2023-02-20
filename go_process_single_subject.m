% ********************************************************************** %
% Title: Preprocessing Script for Go/No-Go Data [Script 2]
% Authors: Armen Bagdasarov & Kenneth Roberts
% Institution: Duke University
% ********************************************************************** %

%% Script Description
% This script runs the preprocessing steps as described below for each
% subjects. This script should not be run, but the file paths need to be
% changed to make sure everything is saved where it should be saved. Script
% 1 (go_loop_over_subjects.m) will call this script.

%% Start
function summary_info = go_process_single_subject(varargin)

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

%% Remove segments without events (i.e., breaks)

% Remove the first 20 seconds of data first (contains SESS, CELL, etc. tags that we do not need)
EEG = pop_select(EEG, 'notime',[0 20]);

% Delete segments of data between event codes if the length of the segment is greater than 4000 ms
% Save 1000 ms of time before the first event code and 1000 ms after the end of the last event code
% This keeps the practice trials, which is fine for now, they will be removed later
EEG  = pop_erplabDeleteTimeSegments(EEG , 'displayEEG',  0, 'endEventcodeBufferMS',  ...
    1000, 'ignoreUseType', 'Ignore', 'startEventcodeBufferMS',  1000, 'timeThresholdMS',  4000 );

%% Filter

% Low-pass filter at 40 Hz and high-pass filter at 1 Hz
EEG = pop_eegfiltnew(EEG, 'hicutoff', 40); 
    % Because low-pass filter is applied at 40 Hz, we do not
    % have to worry too much about potential line noise at 60 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1);

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

% ASR and ICA work best with data filtered at 1 Hz high pass, 
% so let's high pass filter at 1 Hz temporarily
EEG = pop_eegfiltnew(EEG, 'locutoff', 1); 

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

%% ICA
% Extended infomax ICA with PCA dimension reduction
% PCA dimension reduction is necessary because of the large number of 
% channels and relatively short amount of data
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', 50);

% Save ICA plot
ica_plot_path = '[INSERT PATH]'; 
ica_plot_name = [proj.currentId '_go_ica_plot'];
pop_topoplot(EEG, 0, [1:50], 'Independent Components', 0, 'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_plot_path, ica_plot_name), 'png');
close(gcf);

%% Select independent components related to eye or muscle only

% Automatic classification with ICLabel
EEG = pop_iclabel(EEG, 'default');

% Flag components with >= 70% of being eye or muscle
EEG = pop_icflag(EEG, [NaN NaN; 0.7 1; 0.7 1; NaN NaN; NaN NaN; ...
    NaN NaN; NaN NaN]);

% Select components with >= 70% of being eye or muscle
eye_prob = EEG.etc.ic_classification.ICLabel.classifications(:,3);
muscle_prob = EEG.etc.ic_classification.ICLabel.classifications(:,2);
eye_rej = find(eye_prob >= .70);
muscle_rej = find(muscle_prob >= .70);
eye_muscle_rej = [eye_rej; muscle_rej];
eye_muscle_rej = eye_muscle_rej';

% Save retained variance post-ICA in summary info
% This is important because if too much variance is lost, we might want to 
% consider excluding the participant from further analyses.
[projected, pvar] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, EEG.icawinv, eye_muscle_rej);
summary_info.var_retained = 100-pvar;

% Plot only the removed components
ica_rej_plot_path = '[INSERT PATH]'; 
ica_rej_plot_name = [proj.currentId '_go_removed_ics'];
figure % This line is necessary for if there is only 1 component to plot
pop_topoplot(EEG, 0, eye_muscle_rej, 'Independent Components', 0, ...
    'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_rej_plot_path, ica_rej_plot_name), 'png');
close(gcf); close(gcf); % Need to close twice if there is only 1 component to plot

%% Copy EEG ICA fields to EEG_no_rej and remove components with >= 70% of being eye or muscle
% Basically, back-projecting the ICA information from the ASR-reduced data to the full pre-ASR data

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

%% Save file after IC removal but before anything else (i.e., interim results)
interim_path = '[INSERT PATH]'; 
interim_name = [proj.currentId '_go_interim'];
pop_saveset(EEG, fullfile(interim_path, interim_name));

%% Plot channel spectra
% These spectra will likely be noisy because we have not performed
% rejection of epochs yet. Just for visualization of data to note 
% anything wonky.

spectra_time = EEG.xmax * 1000;
figure; pop_spectopo(EEG, 1, [0  spectra_time], 'EEG' , 'freqrange',...
    [2 80],'electrodes','off');

% Save channel spectra plot
spectra_plot_path = '[INSERT PATH]'; 
spectra_plot_name = [proj.currentId '_go_spectra_plot'];
saveas(gca, fullfile(spectra_plot_path, spectra_plot_name), 'png');
close(gcf);

%% Segmentation

% Select only the events that we want; remove all other events like practice trials
EEG = pop_selectevent( EEG, 'type',{'Go_C','Go_I','No_C','No_I','resp'},'deleteevents','on');

% Load event list
el_filename = ['[INSERT PATH]' proj.currentId '_go_event_list.txt'];
EEG  = pop_editeventlist( EEG , 'BoundaryNumeric', { -99}, 'BoundaryString', { 'boundary' }, ...
    'ExportEL', el_filename, 'List', '[INSERT PATH]', ...
    'SendEL2', 'EEG&Text', 'UpdateEEG', 'code', 'Warning', 'off' );

% Upload bins file
el_import_filename =  fullfile('', el_filename);
EEG  = pop_binlister( EEG , 'BDF', '[INSERT PATH]', 'ImportEL', ...
  el_import_filename, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' );

% Extract bin-based epochs and apply baseline correction
EEG = pop_epochbin(EEG , [-500  800],  [-500 -300]); 
% Baseline is -500 to -300 because the ERN can begin prior to the motor response

%% Artifact rejection using TBT plugin

EEG_orig_erp = EEG; % Save the EEG structure

% 1) Max-min amplitude difference
% Basically, the same as ERPLAB's moving window peak-to-peak threshold
% Peak-to-peak amplitude exceeding 100 uV within 200 ms windows sliding by 20 ms
EEG = pop_eegmaxmin(EEG, 1:EEG.nbchan,[-500  796], 100, 200, 20, 0);

% 2) Abnormal values
% Simple voltage thresholding of -150/+150 uV
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -150 , 150 , -0.5 , 0.796, 1, 0); 

% 3) Improbable data
% Based on joint probability, SD = 3 for both local and global thresholds 
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 1, 0, 0);

% Reject based on the epochs selected above
[EEG badlist]= pop_TBT(EEG, EEG.reject.rejmaxminE | EEG.reject.rejthreshE | EEG.reject.rejjpE, 10,1,0); 
    % Do chan/epoch interpolation on ALL 3 types of artifact rejection at once
    % Criteria must be met in at least 10 channels for the epoch to be rejected

% Now for a bit of a hack: We want to copy the interpolated epochs from pop_TBT back into the original structure that we saved
    
    % Use eeg_epochformat to get epoch fields to timelocking event
    [epoch_arr, fn] = eeg_epochformat(EEG.epoch, 'array');
    fn_ind = find(strcmp(fn, 'eventbepoch'));
    eventbepoch_tbt = [epoch_arr{:, fn_ind}];
    
    % Push good events back into the struct
    EEG_orig_erp.data(:,:,eventbepoch_tbt) = EEG.data;
    
    % Set the EEGLAB flags in EEG.reject to reflect epochs flagged by pop_TBT
    EEG_orig_erp.reject.rejmanual = true(1, size(EEG_orig_erp.data, 3));
    EEG_orig_erp.reject.rejmanual(eventbepoch_tbt) = false;
    
    % Sync from EEGLAB to ERPLAB
    % Mark the bad epochs in an ERPLAB-compatible way so ERPLAB knows which epochs are bad
    EEG_orig_erp = pop_syncroartifacts(EEG_orig_erp, 'Direction', 'eeglab2erplab');
    
    % Copy back into EEG
    EEG = EEG_orig_erp;

% Save how many epochs are kept 
summary_info.n_epochs_kept = sum(~EEG.reject.rejmanual);

%% Save final preprocessed files
% EEGLAB-compatible .set format
set_path = '[INSERT PATH]'; 
set_name = [proj.currentId '_go_complete_without_rt_rej']; % This script does NOT include rejection of epochs based on reaction-time data
                                                           % This will require further custom scripts
pop_saveset(EEG, 'filepath', set_path, 'filename', set_name);
