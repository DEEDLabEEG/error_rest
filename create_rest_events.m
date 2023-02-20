% ********************************************************************** %
% Title: Preprocessing Script for Resting-State Data [Script 3]
% Authors: Armen Bagdasarov & Kenneth Roberts
% Institution: Duke University
% ********************************************************************** %

%% Script Description
% This script inserts eyes-open and eyes-closed markers because our data
% does not contain these. This script should not be run because it will be
% caleld in rest_process_single_subject.m.

%% Start
function [EEG, info] = create_rest_events(EEG)

% From rs files with 8 'bgin' events this will:
%   - Relabel those events with rs_open and rs_closed
%   - Generate 'rs_end' events to mark 60s after beginning of resting state
%   - Do checks and warn if any epochs overlap
%   - Do checks and warn if any epoch is too short (goes off end of file)

% Read data out of event
evt_codes = { EEG.event(:).code };
evt_types = { EEG.event(:).type };
evt_lats = [ EEG.event(:).latency ];

% Find 'bgin' events
rs_evt_ind = find(strcmp('bgin', evt_codes));

% Label 'open' and 'closed'
[ EEG.event(rs_evt_ind(1:2:end)).code ] = deal('rs_open');
[ EEG.event(rs_evt_ind(1:2:end)).type ] = deal('rs_open');
[ EEG.event(rs_evt_ind(2:2:end)).code ] = deal('rs_closed');
[ EEG.event(rs_evt_ind(2:2:end)).type ] = deal('rs_closed');

% Add new events that mark the end of each period
nevt = length(EEG.event);
evt_template = EEG.event(rs_evt_ind(1)); % Make a "template" event from first rs event
evt_newlats = [EEG.event(rs_evt_ind).latency] + 60*EEG.srate;

% Ensure no event goes off end of the file, move events to within 10 pnts
% of end of file (leaving a little buffer for resampling/filtering)
evt_newlats_from_end = EEG.pnts - evt_newlats;
evt_newlats = min(evt_newlats, EEG.pnts-10);

[~, si] = sort([evt_lats, evt_newlats]); % Need to sort in case a block is shorter than 60s
EEG.event(si <= nevt) = EEG.event(:); % Spread out old events

% Fill in values of newm 'rs_end' events
[EEG.event(si > nevt)] = deal(EEG.event(rs_evt_ind(1))); % Generate new events by copying first rs 'bgin'
[EEG.event(si > nevt).code] = deal('rs_end'); % Fix events by filling in important fields
[EEG.event(si > nevt).type] = deal('rs_end');
[EEG.event(si > nevt).begintime] = deal('');

evt_newlats = mat2cell(evt_newlats, [1], ones(1, length(evt_newlats)));
[EEG.event(si > nevt).latency] = evt_newlats{:};

% See if blocks overlap
begin_ind = ismember({EEG.event(:).code}, {'rs_open', 'rs_closed'});
end_ind = ismember({EEG.event(:).code}, {'rs_end'});

blocklen = diff([EEG.event(begin_ind).latency EEG.pnts]/EEG.srate);
short_block = find(blocklen < 60);
if ~isempty(short_block)
    warning('Blocks %s too short, length of %s seconds respectively.\n', ...
        mat2str(short_block), mat2str(blocklen(short_block)));
    info.block_overlap = true;
    info.block_int = blocklen;
else
    info.block_overlap = false;
    info.block_int = blocklen;
end

% See if any block has been ended early because it exceeded length of file (or came within 10 points)
if any(evt_newlats_from_end < 10)
    warning('Final %d blocks were truncated by early end of file.', ...
        sum(evt_newlats_from_end <0));
    info.block_truncate = true;
else
    info.block_truncate = false;
end

info.blocklen = ([EEG.event(end_ind).latency] - [EEG.event(begin_ind).latency])/EEG.srate;
