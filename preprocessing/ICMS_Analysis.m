%% ICMS stimulus-triggered averaging analysis

clear; clc; close all;

%% Input files and EMG channels

cerestim = 'CerestimChannelData.mat';
emgFile  = 'EMG_recording.ns4';

% Six tongue muscles included in the analysis
emg_stimulation_channels = [2, 3, 4, 6, 7, 8];

emg_muscles = { ...
    'Right Genioglossus', ...
    'Left Intrinsic', ...
    'Right Intrinsic', ...
    'Left Genioglossus', ...
    'Left Hyoglossus', ...
    'Right Hyoglossus'};

emg_map = dictionary(emg_stimulation_channels, emg_muscles);
tongue_muscles = emg_muscles;

%% 1) Add NPMK path

addpath('path_to_NPMK') % https://github.com/BlackrockNeurotech/NPMK

%% 2) Setup: define number of electrodes, trains, and pulses

numElectrodes       = 96;   % 96-electrode Utah array
numTrainsPerChannel = 34;   % 34 trains per electrode
numPulsesPerTrain   = 15;   % 15 biphasic pulses per train

% For StTA, only the first 14 anodal-phase onsets are used.
% This gives 34 trains x 14 pulses = 476 triggers per electrode.
numAnalyzedPulsesPerTrain = numPulsesPerTrain - 1;

numTrainsTotal = numElectrodes * numTrainsPerChannel;  % 3264 total trains

%% 3) Load the Cerestim channel

CerestimChannelData = load(cerestim).CerestimChannelData;
fs_Utah = 30000;  % 30 kHz
clear NS5

%% 4) Find train start and end indices

threshold = 1e4;
binary_signal = (CerestimChannelData > threshold);

% Rising edges -> train starts
trainStartIndices = find(diff(binary_signal) == 1) + 1;

% Falling edges -> train ends
trainEndIndices = find(diff(binary_signal) == -1) + 1;

if length(trainStartIndices) < numTrainsTotal
    warning('Fewer train starts (%d) than expected (%d)', ...
            length(trainStartIndices), numTrainsTotal);
end

if length(trainEndIndices) < numTrainsTotal
    warning('Fewer train ends (%d) than expected (%d)', ...
            length(trainEndIndices), numTrainsTotal);
end

%% 5) Build structure for all trains

TheTrains = struct();

for t = 1:numTrainsTotal
    channelID = floor((t-1)/numTrainsPerChannel) + 1;  % stimulation order, 1..96
    trainID   = mod((t-1), numTrainsPerChannel) + 1;   % train number, 1..34

    TheTrains(t).ChannelID  = channelID;
    TheTrains(t).TrainID    = trainID;
    TheTrains(t).StartIndex = trainStartIndices(t);
    TheTrains(t).EndIndex   = trainEndIndices(t);
end

% Stimulation order -> physical electrode mapping
% Converts sequential stimulation index to physical electrode number on
% the Utah array grid, based on the electrode map.
stimulation_channels = 1:96;

electrode_channels = [1,  5,  10, 16, 23, 31,  2,  6, ...
                      11, 17, 24, 32,  3,  7, 12, 18, ...
                      25,  4,  8, 13, 19, 26,  9, 14, ...
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

for i = 1:numel(TheTrains)
    TheTrains(i).ElectrodeIDs = channel_map(TheTrains(i).ChannelID);
end

%% 6) Compute pulse timings for each train

% Convert microseconds -> samples at 30 kHz
negSamples   = round(200   * 1e-6 * fs_Utah);  % 200 us cathodal phase
ipiSamples   = round(65535 * 1e-6 * fs_Utah);  % 65.535 ms interphase interval
posSamples   = round(200   * 1e-6 * fs_Utah);  % 200 us anodal phase
pauseSamples = round(53    * 1e-6 * fs_Utah);  % 53 us hardware recovery

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
        pulseEndIdx = posEndIdx   + pauseSamples;

        % Columns:
        % 1 = cathodal phase start
        % 2 = interphase start
        % 3 = anodal phase start
        % 4 = anodal phase end
        % 5 = pulse end
        PulseTimings(p,:) = [negStartIdx, ipiStartIdx, posStartIdx, posEndIdx, pulseEndIdx];

        currentIndex = pulseEndIdx;
    end

    TheTrains(t).PulseTimings = PulseTimings;

    if PulseTimings(end,end) > trainEnd
        overshoot = PulseTimings(end,end) - trainEnd;
        warning('Train %d extends %d samples beyond end edge.', t, overshoot);
    end
end

%% 7) Load EMG data from NS4

openNSx(emgFile);

EMG_Signal_Channels = NS4.Data;
EMG_Signal_Channels = EMG_Signal_Channels(emg_stimulation_channels, :);
fs_EMG = double(NS4.MetaTags.SamplingFreq);  % expected 10,000 Hz

clear NS4

[numEMGChannels, numEMGSamples] = size(EMG_Signal_Channels);

%% 8) Filter the continuous EMG signal and full-wave rectify

% The continuous EMG signal is filtered before stimulus-triggered windows
% are extracted to avoid edge effects from filtering short windows.
hp_cutoff   = 10;
lp_cutoff   = 500;
filterOrder = 4;

nyquistFreq = fs_EMG / 2;
[b, a] = butter(filterOrder, [hp_cutoff/nyquistFreq, lp_cutoff/nyquistFreq]);

pulseSelectionThreshold = 0;  % minimum pulse-level peak z-score for StTA inclusion
stdThreshold = 5;             % final averaged-StTA response classification threshold

Filtered_EMG_Signal_Channels = zeros(size(EMG_Signal_Channels));

for emgCh = 1:numEMGChannels
    Filtered_EMG_Signal_Channels(emgCh,:) = ...
        filtfilt(b, a, double(EMG_Signal_Channels(emgCh,:)));
end

Rectified_EMG_Signal_Channels = abs(Filtered_EMG_Signal_Channels);

%% 9) Extract EMG windows around each analyzed anodal phase

% Aligned to the start of the positive/anodal phase.
% Uses the first 14 anodal-phase onsets per train.
% Window: -20 ms baseline to +40 ms post-trigger analysis.
%
% RawWaveform is extracted from the raw continuous EMG.
% FilteredWave is extracted from the continuously filtered EMG.
% RectifiedWave is extracted from the continuously filtered and rectified EMG.

baselineMs  = 20;  % pre-trigger baseline
analysisMs  = 40;  % post-trigger analysis window
totalMs     = baselineMs + analysisMs;  % 60 ms total

baselineSamples_10k = round(baselineMs * 1e-3 * fs_EMG);
analysisSamples_10k = round(analysisMs * 1e-3 * fs_EMG);
totalSamples_10k    = baselineSamples_10k + analysisSamples_10k;

for t = 1:numTrainsTotal
    PulseTimings = TheTrains(t).PulseTimings;  % 15 x 5
    TheTrains(t).ExtractedEMG = struct([]);

    for p = 1:numAnalyzedPulsesPerTrain  % p = 1..14

        % Align to start of positive/anodal phase.
        referenceIndex_30k = PulseTimings(p, 3);

        % Convert 30 kHz sample index to 10 kHz sample index.
        referenceIndex_10k = round(referenceIndex_30k / 3);

        windowStart_10k = referenceIndex_10k - baselineSamples_10k;
        windowEnd_10k   = referenceIndex_10k + analysisSamples_10k - 1;

        if (windowStart_10k < 1) || (windowEnd_10k > numEMGSamples)
            warning('Train %d, Pulse %d -> out-of-bounds EMG extraction.', t, p);
            continue;
        end

        rawWave_6xN  = EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);
        filtWave_6xN = Filtered_EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);
        rectWave_6xN = Rectified_EMG_Signal_Channels(:, windowStart_10k:windowEnd_10k);

        TheTrains(t).ExtractedEMG(p).RawWaveform  = rawWave_6xN;
        TheTrains(t).ExtractedEMG(p).FilteredWave = filtWave_6xN;
        TheTrains(t).ExtractedEMG(p).RectifiedWave = rectWave_6xN;

        TheTrains(t).ExtractedEMG(p).ReferenceIndex_30k = referenceIndex_30k;
        TheTrains(t).ExtractedEMG(p).ReferenceIndex_10k = referenceIndex_10k;
    end
end

%% 10) Per-electrode stimulus-triggered averaging

% Each individual rectified EMG segment is evaluated using its own
% 20 ms pre-trigger baseline and 40 ms post-trigger peak.
% A segment is included only when its pulse-level peak z-score is greater
% than or equal to pulseSelectionThreshold.
% Only selected rectified EMG segments are averaged.
%
% Final StTA response magnitude:
%
%   sdAbove = (peak post-trigger EMG - baseline mean) / baseline SD
%
% The final peak is taken from the 40 ms post-trigger window only.
% The final baseline mean and SD are taken from the 20 ms pre-trigger
% portion of the averaged selected waveform.

timeZeroSample = baselineSamples_10k + 1;

analysisStartSample = timeZeroSample;
analysisEndSample   = baselineSamples_10k + analysisSamples_10k;

ElectrodeStTA = struct();

for e = 1:numElectrodes

    % Match by physical electrode number, not stimulation-order number.
    theseTrains = find([TheTrains.ElectrodeIDs] == e);

    validRectified = cell(numEMGChannels, 1);

    for idxT = 1:length(theseTrains)

        t = theseTrains(idxT);
        numPulses = length(TheTrains(t).ExtractedEMG);

        for p = 1:numPulses

            if ~isfield(TheTrains(t).ExtractedEMG(p), 'RectifiedWave')
                continue;
            end

            rectWave_6xN = TheTrains(t).ExtractedEMG(p).RectifiedWave;

            for emgCh = 1:numEMGChannels

                waveCh = rectWave_6xN(emgCh,:);

                % Evaluate each pulse waveform using its own pre-trigger
                % baseline and post-trigger peak.
                pulseBaselineWave = waveCh(1:baselineSamples_10k);
                pulseBaselineMean = mean(pulseBaselineWave);
                pulseBaselineSD   = std(pulseBaselineWave);

                pulseAnalysisWave = waveCh(analysisStartSample:analysisEndSample);
                pulsePeakAmplitude = max(pulseAnalysisWave);

                % A zero baseline SD cannot produce a valid pulse z-score.
                if pulseBaselineSD <= 0
                    continue;
                end

                pulseSdAbove = ...
                    (pulsePeakAmplitude - pulseBaselineMean) / pulseBaselineSD;

                % Include only pulse waveforms meeting the selection threshold.
                if pulseSdAbove < pulseSelectionThreshold
                    continue;
                end

                if isempty(validRectified{emgCh})
                    validRectified{emgCh} = waveCh(:);
                else
                    validRectified{emgCh} = [validRectified{emgCh}, waveCh(:)];
                end
            end
        end
    end

    ElectrodeStTA(e).ElectrodeID = e;
    ElectrodeStTA(e).Channel = struct();

    for emgCh = 1:numEMGChannels

        ElectrodeStTA(e).Channel(emgCh).MuscleName = emg_muscles{emgCh};

        if isempty(validRectified{emgCh})

            ElectrodeStTA(e).Channel(emgCh).MeanRectified = [];
            ElectrodeStTA(e).Channel(emgCh).Results       = [];

            continue;
        end

        % Rows = time samples; columns = valid triggers.
        waveMatrix = validRectified{emgCh};

        % StTA waveform: average threshold-selected rectified segments.
        meanRectWave = mean(waveMatrix, 2)';   % 1 x totalSamples_10k

        ElectrodeStTA(e).Channel(emgCh).MeanRectified = meanRectWave;

        % Baseline from 20 ms pre-trigger portion of averaged waveform.
        baselineWave = meanRectWave(1:baselineSamples_10k);
        muBase = mean(baselineWave);
        sdBase = std(baselineWave);

        % Peak from 40 ms post-trigger analysis window only.
        analysisWave = meanRectWave(analysisStartSample:analysisEndSample);
        [peakAmplitude, peakIdxRel] = max(analysisWave);
        peakIdx = analysisStartSample + peakIdxRel - 1;

        if sdBase > 0
            sdAbove = (peakAmplitude - muBase) / sdBase;
        else
            sdAbove = NaN;
            warning('Baseline SD is zero for electrode %d, EMG channel %d.', e, emgCh);
        end

        % Optional response classification.
        % This does not affect the StTA waveform or sdAbove value.
        effectClass = 'none';

        if ~isnan(sdAbove) && sdAbove >= stdThreshold
            effectClass = 'mod/strong';
        end

        % Latency relative to the anodal trigger.
        peakLatency_ms = ((peakIdx - timeZeroSample) / fs_EMG) * 1000;

        results.BaselineMean       = muBase;
        results.BaselineSD         = sdBase;
        results.PeakSample         = peakIdx;
        results.PeakLatency_ms     = peakLatency_ms;
        results.PeakAmplitude      = peakAmplitude;
        results.EffectClass        = effectClass;
        results.sdAbove            = sdAbove;
        results.NumTriggers        = size(waveMatrix, 2);
        results.PulseSelectionThreshold = pulseSelectionThreshold;

        % Onset/offset are not estimated because response magnitude is
        % defined by the post-trigger peak z-score.
        results.OnsetSample       = NaN;
        results.OffsetSample      = NaN;
        results.OnsetLatency_ms   = NaN;
        results.OffsetLatency_ms  = NaN;

        ElectrodeStTA(e).Channel(emgCh).Results = results;
    end
end

%% 11) Grid StTA: SD above baseline per physical electrode per muscle

gridSTA = struct();

for i = 1:numElectrodes

    gridSTA(i).ElectrodeID = i;

    for muscle = 1:numEMGChannels

        rawName    = string(emg_muscles(muscle));
        muscleName = matlab.lang.makeValidName(rawName);

        muscleResult = ElectrodeStTA(i).Channel(muscle).Results;

        if isempty(muscleResult) || ~isfield(muscleResult, 'sdAbove')
            gridSTA(i).(muscleName) = NaN;
        else
            gridSTA(i).(muscleName) = muscleResult.sdAbove;
        end
    end
end

%% 12) Binary response maps based on optional 5-SD classification

% staMap(i,muscle) = 1 if the StTA peak is >= 5 baseline SD.
% This is only a thresholded summary map. The continuous map value is
% stored in gridSTA as sdAbove.

staMap = zeros(numElectrodes, numEMGChannels);

for i = 1:numElectrodes

    for muscle = 1:numEMGChannels

        muscleResult = ElectrodeStTA(i).Channel(muscle).Results;

        if isempty(muscleResult) || ~isfield(muscleResult, 'EffectClass')
            staMap(i, muscle) = 0;
        elseif strcmp(muscleResult.EffectClass, 'mod/strong')
            staMap(i, muscle) = 1;
        else
            staMap(i, muscle) = 0;
        end
    end
end

%% Optional save

% save('StTA_Output.mat', ...
%     'TheTrains', ...
%     'ElectrodeStTA', ...
%     'gridSTA', ...
%     'staMap', ...
%     'emg_muscles', ...
%     'tongue_muscles', ...
%     'electrode_channels', ...
%     'baselineMs', ...
%     'analysisMs', ...
%     'fs_EMG', ...
%     'fs_Utah');
