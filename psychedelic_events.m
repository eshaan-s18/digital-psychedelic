clear;
eeglab;
close;
clc

% We have three streams per recording
% EEG, event markets for state transitions, bci output from neuropype
% (alpha) for state transitions

%% Loading Dataset
folder = 'digital_psychedelic';
fp = fullfile('data', '2305_dmn_1372_condition.xdf');

[~, name, ext] = fileparts(fp);
filename = name;
output_folder = 'sets/';

streams = load_xdf(fp);

for i = 1:length(streams)
    s = streams{i};
    fprintf('steam %i (%s): $s\n', i, s.info.type, s.info.name)
end

%% Extracting Streams

streamState = streams{cellfun(@(x) strcmp(x.info.name, 'StateMarkers'),streams)};
streamEEG = streams{cellfun(@(x) strcmp(x.info.type, 'EEG'),streams)};

%% Loading EEG in eeglab

EEG = pop_loadxdf(fp);

%% Inserting Events
event_types = streamState.time_series;
event_latency = num2cell((streamState.time_stamps-streamEEG.time_stamps(1))*EEG.srate);

[EEG.event(1:length(event_types)).type] = event_types{:};
[EEG.event(1:length(event_latency)).latency] = event_latency{:};

% -----------------------------------------------------------------------
%% Preprocess

eeglab redraw

%{

This pipeline is as followed:

1. Load in XDF
2. Load Channel Information using an external file
3. High Pass Filter at 1Hz
4. Low Pass Filter at 128Hz
5. Notch Filter at 60Hz
6. Resample to 256Hz
7. Remove bad channels and apply ASR
    a. Do not use burst reject and use ASR reconstruction
    b. Save what channels were deleted as a field named bad_channels
8. Interpolate bad channels
    a. Insert a zero channel (as we have discussed) to perform average
    re-reference

%}

%% 2. Load Channel Info.
nonEEG_labels = {
    'Trig1', 'EX1', 'EX2', 'EX3', 'EX4', 'EX5', 'EX6', 'EX7', 'EX8'
    };

% Removing nonEEG chans
EEG = pop_select(EEG, 'nochannel', nonEEG_labels);

% Loading chan info
hp = fullfile('headset_channel_locations', 'bio_semi_64.ced');
EEG = pop_chanedit(EEG, 'load', hp);


%% 3. High Pass Filter at 1Hz
EEG = pop_eegfiltnew( ...
    EEG, ...
    'locutoff', 1, ...
    'hicutoff', 0);


%% 4. Low Pass Filter
EEG = pop_eegfiltnew( ...
    EEG, ...
    'locutoff', 0, ...
    'hicutoff', 128);


%% 5. Notch Filter at 60Hz
EEG = pop_eegfiltnew( ...
    EEG, ...
    'locutoff', 59, ...
    'hicutoff', 61, ...
    'revfilt', 1);


%% 6. Resample data at 256Hz
EEG = pop_resample(EEG, 256);


%% 7. Remove bad channels and apply ASR
    
    % Save original EEG data
    originalEEG = EEG;
    
    % Clean Raw Data
    EEG = pop_clean_rawdata( ...
        originalEEG, ...
        'FlatlineCriterion', 5, ...
        'ChannelCriterion', 0.8, ...
        'LineNoiseCriterion', 4, ...
        'Highpass', 'off', ...
        'BurstCriterion', 20, ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off', ...
        'Distance', 'Euclidian');
    
    % Creating a field in the EEG struct to store bad channels
    new_labels = {EEG.chanlocs.labels};
    org_labels = {originalEEG.chanlocs.labels};
    [~,ia,~] = intersect(org_labels, new_labels, 'stable');
    removed_channels = org_labels;
    removed_channels(ia) = [];
    EEG.bad_channels = removed_channels;

%% 9. Interpolating
    
% spherical spline interpolation
EEG = pop_interp(EEG, originalEEG.chanlocs);


%% 10. Re-reference

reference_chan = {};

if not(isempty(reference_chan))
    EEG = pop_reref(EEG, reference_chan);
else
    EEG.nbchan = EEG.nbchan + 1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1, EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select(EEG, 'nochannel', {'initialReference'});
end

% -----------------------------------------------------------------------

%% Exploring events
idx_st = find(contains({EEG.event.type},{'StateTransition'}));
% Event Information
st_lat = [EEG.event(idx_st).latency];
st_evs = {EEG.event(idx_st).type};
st_dur = diff([EEG.event(idx_st).latency])/EEG.srate;


%% Epoching
freqs = 0:128;
% boolean array
ba = (freqs >= 8) & (freqs <= 12); % !review these values!

% Initialize empty structure to hold epochs
epoch_st = cell(length(st_evs)-1,1);
for i = 1:length(st_evs)-1
    ev = st_evs{i};
    dur = st_dur(i);
    epoch_st{i} = pop_epoch(EEG, {ev}, [0,dur]);

    epoch_name = sprintf('%s_epoch_%d.set', filename, i);
    pop_saveset(epoch_st{i}, 'filename', epoch_name, 'filepath', output_folder);

end

%% Save dataset
save(fullfile(output_folder, '2305_dmn_1372_condition_epoch_st.mat'), 'epoch_st');




