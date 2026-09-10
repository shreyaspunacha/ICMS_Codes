%% M1 Utah Array ICMS Stimulation
% Requires the CereStim MATLAB API and compatible stimulation hardware.

clear;
close all;
clc;

%% Stimulation parameters
% Balanced biphasic ICMS pulse
% 0.2 ms cathodal phase followed by 0.2 ms anodal phase
% Low-intensity stimulation (<20 uA)
% 15 pulses per train
% 34 trains per electrode (510 pulses total)
% Inter-electrode interval: 0.5 s

%% Add path to the CereStim API
addpath(genpath('path/to/CereStim-API'));

%% Define stimulator object
stimulator = cerestim96();

%% Find and connect to an available stimulator
Serials = stimulator.scanForDevices();

if ~isempty(Serials)
    stimulator.selectDevice(Serials(1));
    stimulator.connect();
end

%% Define stimulation parameters
stimulator.setStimPattern( ...
    'waveform', 1, ...
    'polarity', 1, ...
    'pulses', 15, ...
    'amp1', 18, ...
    'amp2', 18, ...
    'width1', 200, ...
    'width2', 200, ...
    'interphase', 65535, ...
    'frequency', 15);

%% Stimulate sequentially across the three electrode banks

channels = 1:32;

for i_bank = 1:3

    if i_bank == 1

        A_chans = channels;
        fprintf('Moving to Bank A\n');

        for i = 1:length(A_chans)

            ChannelNumber = A_chans(i);

            stimulator.beginSequence;
            stimulator.autoStim(ChannelNumber, 1);
            pause(0.5);
            stimulator.endSequence;

            fprintf('Stimulating Channel Number = %d\n', ChannelNumber);

            for iteration = 1:34
                stimulator.play(iteration);
                fprintf('Iteration = %d\n', iteration);
                stimulator.stop();
            end
        end

    elseif i_bank == 2

        B_chans = 32 + channels;
        fprintf('Moving to Bank B\n');

        for i = 1:length(B_chans)

            ChannelNumber = B_chans(i);

            stimulator.beginSequence;
            stimulator.autoStim(ChannelNumber, 1);
            pause(0.5);
            stimulator.endSequence;

            fprintf('Stimulating Channel Number = %d\n', ChannelNumber);

            for iteration = 1:34
                stimulator.play(iteration);
                fprintf('Iteration = %d\n', iteration);
                stimulator.stop();
            end
        end

    else

        C_chans = 64 + channels;
        fprintf('Moving to Bank C\n');

        for i = 1:length(C_chans)

            ChannelNumber = C_chans(i);

            stimulator.beginSequence;
            stimulator.autoStim(ChannelNumber, 1);
            pause(0.5);
            stimulator.endSequence;

            fprintf('Stimulating Channel Number = %d\n', ChannelNumber);

            for iteration = 1:34
                stimulator.play(iteration);
                fprintf('Iteration = %d\n', iteration);
                stimulator.stop();
            end
        end
    end
end

%% Save stimulation parameters
save('M1Params.mat');

%% Disconnect the stimulator
stimulator.disconnect();
