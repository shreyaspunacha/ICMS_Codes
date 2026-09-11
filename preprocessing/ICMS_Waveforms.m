%% ========================================================================
% ICMS StTA PER-STIMULATION WAVEFORM SELECTION
%
% Purpose:
%   Evaluate each individual rectified pulse waveform using its own
%   pre-trigger baseline and retain waveforms whose analysis-window peak
%   meets the pulse-selection threshold (set to zero here) defined by 
%   stdThreshold.
%
% Main output for subsequent waveform-level statistics:
%   DRMS_Input.SelectedWaveforms{electrode, muscle}
%
% Each cell contains:
%   time samples x number of threshold-passing pulse waveforms
%
% Threshold condition:
%   analysisPeak >= baselineMean + stdThreshold * baselineSD
%
% In this script, stdThreshold = 0.
%
% Notes:
%   - 34 trains per electrode
%   - 15 pulses per train
%   - first 14 pulses are analyzed
%   - 34 x 14 = 476 candidate stimulations/electrode
%
% Run once per condition:
%   Control
%   AllNB
%   SelectiveNB
%
%% ========================================================================

clear; clc; close all;
tic;

%% -------------------- USER SETTINGS -------------------------------------

% Condition label
conditionLabel = 'Control';   % Change to 'Control', 'AllNB', or 'SelectiveNB'

% Input files
cerestim = 'CerestimChannelData.mat';
emgFile  = 'EMG_recording.ns4';

% Output file
outputFile = ['DRMS_Input_' conditionLabel '_SelectedWaveforms.mat'];

% NPMK path
% addpath('path/to/NPMK');

% Number settings
numElectrodes       = 96;
numTrainsPerChannel = 34;
numPulsesPerTrain   = 15;
pulsesPerTrainAnalyzed = numPulsesPerTrain - 1;  % 14
numStimPerElectrode = numTrainsPerChannel * pulsesPerTrainAnalyzed; % 476
numTrainsTotal      = numElectrodes * numTrainsPerChannel;

% EMG channel/muscle names
% Channel indices correspond to the six tongue-muscle channels used here.
emg_stimulation_channels = [2, 3, 4, 6, 7, 8];

emg_muscles = { ...
    'Right Genioglossus', ...
    'Left Intrinsic', ...
    'Right Intrinsic', ...
    'Left Genioglossus', ...
    'Left Hyoglossus', ...
    'Right Hyoglossus'};

numMuscles = numel(emg_muscles);

%% -------------------- LOAD CERESTIM CHANNEL -----------------------------

CerestimChannelData = load(cerestim).CerestimChannelData;
fs_Utah = 30000;

%% -------------------- FIND TRAIN START/END INDICES ----------------------

threshold = 1e4;
binary_signal = CerestimChannelData > threshold;

trainStartIndices = find(diff(binary_signal) == 1) + 1;
trainEndIndices   = find(diff(binary_signal) == -1) + 1;

if length(trainStartIndices) < numTrainsTotal
    warning('Fewer train starts (%d) than expected (%d).', ...
        length(trainStartIndices), numTrainsTotal);
end

if length(trainEndIndices) < numTrainsTotal
    warning('Fewer train ends (%d) than expected (%d).', ...
        length(trainEndIndices), numTrainsTotal);
end

%% -------------------- BUILD TRAIN STRUCTURE -----------------------------

TheTrains = struct();

for t = 1:numTrainsTotal

    channelID = floor((t-1) / numTrainsPerChannel) + 1;  % 1 to 96
    trainID   = mod((t-1), numTrainsPerChannel) + 1;     % 1 to 34

    TheTrains(t).ChannelID  = channelID;
    TheTrains(t).TrainID    = trainID;
    TheTrains(t).StartIndex = trainStartIndices(t);
    TheTrains(t).EndIndex   = trainEndIndices(t);
end

%% -------------------- STIMULATION CHANNEL TO ELECTRODE ID MAP -----------

stimulation_channels = 1:96;

electrode_channels = [ ...
    1, 5, 10, 16, 23, 31, 2, 6, ...
    11, 17, 24, 32, 3, 7, 12, 18, ...
    25, 4, 8, 13, 19, 26, 9, 14, ...
    20, 27, 15, 21, 28, 22, 29, 30, ...
    58, 49, 40, 59, 50, 41, 60, 51, ...
    42, 33, 61, 52, 43, 34, 62, 53, ...
    44, 35, 63, 54, 45, 36, 64, 55, ...
    46, 37, 56, 47, 38, 57, 48, 39, ...
    67, 76, 85, 94, 68, 77, 86, 95, ...
    69, 78, 87, 96, 70, 79, 88, 71, ...
    80, 89, 72, 81, 90, 73, 82, 91, ...
    65, 74, 83, 92, 66, 75, 84, 93];

channel_map = dictionary(stimulation_channels, electrode_channels);

for t = 1:numel(TheTrains)
    TheTrains(t).ElectrodeID = channel_map(TheTrains(t).ChannelID);
end

%% -------------------- COMPUTE PULSE TIMINGS -----------------------------

negSamples   = round(200   * 1e-6 * fs_Utah);
ipiSamples   = round(65535 * 1e-6 * fs_Utah);
posSamples   = round(200   * 1e-6 * fs_Utah);
pauseSamples = round(53    * 1e-6 * fs_Utah);

for t = 1:numTrainsTotal

    trainStart = TheTrains(t).StartIndex;
    trainEnd   = TheTrains(t).EndIndex;

    PulseTimings = zeros(numPulsesPerTrain, 5);
    currentIndex = trainStart;

    for p = 1:numPulsesPerTrain

        negStartIdx = currentIndex;
        ipiStartIdx = negStartIdx + negSamples;
        posStartIdx = ipiStartIdx + ipiSamples;
        posEndIdx   = posStartIdx + posSamples;
        pulseEndIdx = posEndIdx + pauseSamples;

        PulseTimings(p,:) = [ ...
            negStartIdx, ...
            ipiStartIdx, ...
            posStartIdx, ...
            posEndIdx, ...
            pulseEndIdx];

        currentIndex = pulseEndIdx;
    end

    TheTrains(t).PulseTimings = PulseTimings;

    if PulseTimings(end,end) > trainEnd
        overshoot = PulseTimings(end,end) - trainEnd;
        warning('Train %d extends %d samples beyond end edge.', t, overshoot);
    end
end

%% -------------------- LOAD EMG DATA -------------------------------------

openNSx(emgFile);

EMG_Signal_Channels = NS4.Data;
EMG_Signal_Channels = EMG_Signal_Channels(emg_stimulation_channels, :);
fs_EMG = double(NS4.MetaTags.SamplingFreq);

clear NS4

[numEMGChannels, numEMGSamples] = size(EMG_Signal_Channels);

if numEMGChannels ~= numMuscles
    warning('Number of EMG channels (%d) does not match number of muscle labels (%d).', ...
        numEMGChannels, numMuscles);
end

%% -------------------- FILTER CONTINUOUS EMG AND RECTIFY -----------------

hp_cutoff   = 10;
lp_cutoff   = 500;
filterOrder = 4;

nyquistFreq = fs_EMG / 2;
[b,a] = butter(filterOrder, [hp_cutoff/nyquistFreq, lp_cutoff/nyquistFreq]);

Filtered_EMG_Signal_Channels = zeros(size(EMG_Signal_Channels));

for emgCh = 1:numEMGChannels
    Filtered_EMG_Signal_Channels(emgCh,:) = ...
        filtfilt(b, a, double(EMG_Signal_Channels(emgCh,:)));
end

Rectified_EMG_Signal_Channels = abs(Filtered_EMG_Signal_Channels);

%% -------------------- EXTRACT EMG WINDOWS -------------------------------

baselineMs = 20;
analysisMs = 40;

baselineSamples_10k = round(baselineMs * 1e-3 * fs_EMG);
analysisSamples_10k = round(analysisMs * 1e-3 * fs_EMG);
totalSamples_10k    = baselineSamples_10k + analysisSamples_10k;

for t = 1:numTrainsTotal

    PulseTimings = TheTrains(t).PulseTimings;
    TheTrains(t).ExtractedEMG = struct([]);

    for p = 1:pulsesPerTrainAnalyzed

        referenceIndex_30k = PulseTimings(p, 3);  % positive phase start
        referenceIndex_10k = round(referenceIndex_30k / 3);

        windowStart_10k = referenceIndex_10k - baselineSamples_10k;
        windowEnd_10k   = referenceIndex_10k + analysisSamples_10k - 1;

        if windowStart_10k < 1 || windowEnd_10k > numEMGSamples
            warning('Train %d, pulse %d out of bounds.', t, p);
            continue;
        end

        rawWave  = EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);
        filtWave = Filtered_EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);
        rectWave = Rectified_EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);

        TheTrains(t).ExtractedEMG(p).RawWaveform      = rawWave;
        TheTrains(t).ExtractedEMG(p).FilteredWave     = filtWave;
        TheTrains(t).ExtractedEMG(p).RectifiedWave    = rectWave;
    end
end

%% -------------------- APPLY PULSE-SELECTION THRESHOLD --------------------

stdThreshold = 0;  % pulse-level peak threshold in baseline SD

% Per-pulse values indexed by actual mapped electrode ID.
PulseSdAbovePeak = nan(numElectrodes, numStimPerElectrode, numMuscles);
PulseSdAboveMean = nan(numElectrodes, numStimPerElectrode, numMuscles);
PulseResponse    = nan(numElectrodes, numStimPerElectrode, numMuscles);

% Full threshold-passing rectified waveforms for subsequent statistics.
% Each cell is: time samples x number of selected pulse waveforms.
SelectedWaveforms = cell(numElectrodes, numMuscles);

% Pulse identity and pulse-level statistics for each retained waveform.
SelectedPulseInfo = cell(numElectrodes, numMuscles);

% Optional: same per-pulse values indexed by stimulation channel order.
PulseSdAbovePeak_ChannelID = nan(numElectrodes, numStimPerElectrode, numMuscles);
PulseSdAboveMean_ChannelID = nan(numElectrodes, numStimPerElectrode, numMuscles);
PulseResponse_ChannelID    = nan(numElectrodes, numStimPerElectrode, numMuscles);

for t = 1:numTrainsTotal

    channelID  = TheTrains(t).ChannelID;
    electrodeID = TheTrains(t).ElectrodeID;
    trainID    = TheTrains(t).TrainID;

    numPulses = length(TheTrains(t).ExtractedEMG);

    for p = 1:numPulses

        if ~isfield(TheTrains(t).ExtractedEMG(p), 'RawWaveform')
            continue;
        end

        rectWave = double(TheTrains(t).ExtractedEMG(p).RectifiedWave);
        numCh = size(rectWave, 1);

        stimIdx = (trainID - 1) * pulsesPerTrainAnalyzed + p;

        for emgCh = 1:numCh

            baselineWave = rectWave(emgCh, 1:baselineSamples_10k);

            % analysis starts after baseline
            analysisWave = rectWave(emgCh, baselineSamples_10k+1:end);

            baseMean = mean(baselineWave);
            sdBase   = std(baselineWave);

            analysisMean = mean(analysisWave);
            analysisPeak = max(analysisWave);

            responded = sdBase > 0 && analysisPeak >= (baseMean + stdThreshold * sdBase);

            if sdBase > 0
                sdAboveMean = (analysisMean - baseMean) / sdBase;
                sdAbovePeak = (analysisPeak - baseMean) / sdBase;
            else
                sdAboveMean = NaN;
                sdAbovePeak = NaN;
            end

            % Store in TheTrains
            TheTrains(t).ExtractedEMG(p).Response(emgCh) = responded;

            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).BaselineMean = baseMean;
            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).BaselineSD   = sdBase;
            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).AnalysisMean = analysisMean;
            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).AnalysisPeak = analysisPeak;
            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).SdAboveMean  = sdAboveMean;
            TheTrains(t).ExtractedEMG(p).PulseStats(emgCh).SdAbovePeak  = sdAbovePeak;

            % Store indexed by mapped electrode ID
            PulseSdAbovePeak(electrodeID, stimIdx, emgCh) = sdAbovePeak;
            PulseSdAboveMean(electrodeID, stimIdx, emgCh) = sdAboveMean;
            PulseResponse(electrodeID, stimIdx, emgCh)    = double(responded);

            % Store indexed by stimulation channel ID
            PulseSdAbovePeak_ChannelID(channelID, stimIdx, emgCh) = sdAbovePeak;
            PulseSdAboveMean_ChannelID(channelID, stimIdx, emgCh) = sdAboveMean;
            PulseResponse_ChannelID(channelID, stimIdx, emgCh)    = double(responded);

            % Retain only waveforms that pass the peak-based selection threshold.
            if responded
                selectedWave = rectWave(emgCh,:).';

                if isempty(SelectedWaveforms{electrodeID, emgCh})
                    SelectedWaveforms{electrodeID, emgCh} = selectedWave;
                else
                    SelectedWaveforms{electrodeID, emgCh}(:,end+1) = selectedWave;
                end

                selectedIndex = size(SelectedWaveforms{electrodeID, emgCh}, 2);

                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).ChannelID    = channelID;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).ElectrodeID  = electrodeID;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).TrainID      = trainID;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).PulseID      = p;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).StimIdx      = stimIdx;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).BaselineMean = baseMean;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).BaselineSD   = sdBase;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).AnalysisPeak = analysisPeak;
                SelectedPulseInfo{electrodeID, emgCh}(selectedIndex).SdAbovePeak  = sdAbovePeak;
            end
        end
    end
end

% Number retained and fraction relative to the nominal candidate pulses per electrode-muscle pair.
NumSelected = zeros(numElectrodes, numMuscles);

for electrodeID = 1:numElectrodes
    for emgCh = 1:numMuscles
        NumSelected(electrodeID, emgCh) = ...
            size(SelectedWaveforms{electrodeID, emgCh}, 2);
    end
end

SelectionFraction = NumSelected / numStimPerElectrode;

%% -------------------- CREATE d_rms INPUT STRUCT -------------------------

DRMS_Input = struct();

DRMS_Input.ConditionLabel = conditionLabel;

% Main waveform-level input for subsequent statistics.
% Each cell contains: time samples x selected pulse waveforms.
DRMS_Input.Data = SelectedWaveforms;

DRMS_Input.DataType = 'ThresholdSelectedRectifiedWaveforms';
DRMS_Input.ElectrodeIndexMode = 'MappedElectrodeID';
DRMS_Input.Dimensions = 'cell array: electrode x muscle; each cell = time x selected pulses';

DRMS_Input.SelectedWaveforms = SelectedWaveforms;
DRMS_Input.SelectedPulseInfo = SelectedPulseInfo;
DRMS_Input.NumSelected = NumSelected;
DRMS_Input.SelectionFraction = SelectionFraction;
DRMS_Input.PassMask = (PulseResponse == 1);

DRMS_Input.PulseSdAbovePeak = PulseSdAbovePeak;
DRMS_Input.PulseSdAboveMean = PulseSdAboveMean;
DRMS_Input.PulseResponse    = PulseResponse;

DRMS_Input.PulseSdAbovePeak_ChannelID = PulseSdAbovePeak_ChannelID;
DRMS_Input.PulseSdAboveMean_ChannelID = PulseSdAboveMean_ChannelID;
DRMS_Input.PulseResponse_ChannelID    = PulseResponse_ChannelID;

DRMS_Input.emg_muscles = emg_muscles;
DRMS_Input.numElectrodes = numElectrodes;
DRMS_Input.numMuscles = numMuscles;
DRMS_Input.numTrainsPerChannel = numTrainsPerChannel;
DRMS_Input.numPulsesPerTrain = numPulsesPerTrain;
DRMS_Input.pulsesPerTrainAnalyzed = pulsesPerTrainAnalyzed;
DRMS_Input.numStimPerElectrode = numStimPerElectrode;
DRMS_Input.fs_EMG = fs_EMG;
DRMS_Input.fs_Utah = fs_Utah;
DRMS_Input.baselineMs = baselineMs;
DRMS_Input.analysisMs = analysisMs;
DRMS_Input.baselineSamples_10k = baselineSamples_10k;
DRMS_Input.analysisSamples_10k = analysisSamples_10k;
DRMS_Input.stdThreshold = stdThreshold;
DRMS_Input.channel_map = channel_map;

%% -------------------- QUICK SANITY CHECKS -------------------------------

fprintf('\nFinished condition: %s\n', conditionLabel);
fprintf('DRMS_Input.Data size: %d x %d cell array\n', size(DRMS_Input.Data));
fprintf('Expected size: 96 x 6 cells\n');

nSelectedTotal = sum(NumSelected(:));
fprintf('Threshold-passing waveforms retained: %d\n', nSelectedTotal);
fprintf('Candidate waveforms: %d\n', ...
    numElectrodes * numStimPerElectrode * numMuscles);

%% -------------------- SAVE OUTPUT ---------------------------------------

save(outputFile, ...
    'DRMS_Input', ...
    'SelectedWaveforms', ...
    'SelectedPulseInfo', ...
    'NumSelected', ...
    'SelectionFraction', ...
    'PulseSdAbovePeak', ...
    'PulseSdAboveMean', ...
    'PulseResponse', ...
    'PulseSdAbovePeak_ChannelID', ...
    'PulseSdAboveMean_ChannelID', ...
    'PulseResponse_ChannelID', ...
    'TheTrains', ...
    'emg_muscles', ...
    'conditionLabel', ...
    'numElectrodes', ...
    'numTrainsPerChannel', ...
    'pulsesPerTrainAnalyzed', ...
    'numStimPerElectrode', ...
    '-v7.3');

fprintf('Saved output file:\n%s\n', outputFile);

toc;