%% ========================================================================
% Plot StTA heatmaps across three conditions
%
% One figure is generated per muscle with three subplots:
%   Control | Combined (All) NB | Selective NB
%
% A common colorbar maximum is used across the three conditions for each
% muscle so that response magnitudes are visually comparable.
%
% Each MAT file must contain gridSTA, a 1 x 96 structure containing
% ElectrodeID and the corresponding muscle fields.
%
% Each heatmap cell represents StTA response magnitude in SD above baseline.
%% ========================================================================

clear; clc; close all;

%% -------------------- USER SETTINGS -------------------------------------

% Subject/session label
animalLabel = 'Subject';

% Condition files and labels
conditionFiles = { ...
    'Subject_Control_GridSTA.mat', ...
    'Subject_AllNB_GridSTA.mat', ...
    'Subject_SelectiveNB_GridSTA.mat'};

conditionLabels = { ...
    'Control', ...
    'All NB', ...
    'Selective NB'};

numConditions = numel(conditionFiles);

% Save figures?
saveFigures = true;
saveFolder = fullfile(pwd, 'StTA_Heatmaps_ThreeConditionComparison');

if saveFigures && ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Plot only the six tongue muscles included in the analysis
plotAllMuscles = false;

% Muscle fields must match the field names in gridSTA.
musclesToPlot = { ...
    'RightHyoglossus', ...
    'LeftIntrinsic', ...
    'RightIntrinsic', ...
    'RightGenioglossus', ...
    'LeftHyoglossus', ...
    'LeftGenioglossus'};

% -------------------------------------------------------------------------
% COLOR SCALE SETTINGS
% -------------------------------------------------------------------------
% Available setting:
%   trueMaxPerMuscle
%
% This uses the actual maximum SD-above-baseline value across the three
% conditions for each muscle.
%
% Alternative:
%   robust99PerMuscle
%
% This uses the 99th percentile across the three conditions for each muscle.
% Use this only if one extreme electrode compresses the rest of the map.
colorScaleMode = 'trueMaxPerMuscle';

% Lower colorbar limit
lowerCLim = 0;

% Minimum upper limit in case all values are very small or invalid
minimumUpperCLim = 1;

% Figure style
fontName = 'Times New Roman';
figurePosition = [100 100 1800 650];

mainTitleFontSize = 26;
subplotTitleFontSize = 24;
numberFontSize = 13;
colorbarFontSize = 22;
labelFontSize = 24;

% Arrow settings
showArrow = true;
arrowColor = [0 0 1];   % blue
arrowLineWidth = 3;

% Save format
savePNG = true;
savePDF = true;

%% -------------------- LOAD ALL THREE CONDITIONS -------------------------

gridSTA_all = cell(1, numConditions);

for k = 1:numConditions

    if ~isfile(conditionFiles{k})
        error('File not found: %s', conditionFiles{k});
    end

    S = load(conditionFiles{k});

    if ~isfield(S, 'gridSTA')
        error('File %s does not contain variable gridSTA.', conditionFiles{k});
    end

    gridSTA_all{k} = S.gridSTA;

    if numel(gridSTA_all{k}) ~= 96
        warning('File %s contains %d electrodes instead of 96.', ...
            conditionFiles{k}, numel(gridSTA_all{k}));
    end
end

%% -------------------- GET MUSCLE FIELDS ---------------------------------

if plotAllMuscles

    % Use fields from the first condition as reference
    allFields = fieldnames(gridSTA_all{1});

    % Remove non-muscle or bookkeeping fields
    nonMuscleFields = { ...
        'ElectrodeID', ...
        'ChannelID', ...
        'TrainID', ...
        'StartIndex', ...
        'EndIndex'};

    candidateFields = setdiff(allFields, nonMuscleFields, 'stable');

    % Keep only fields that exist in all conditions and contain numeric data
    musclesToPlot = {};

    for f = 1:numel(candidateFields)

        fieldName = candidateFields{f};
        validField = true;

        for k = 1:numConditions
            if ~isfield(gridSTA_all{k}, fieldName)
                validField = false;
                break;
            end

            vals = getGridSTAValues(gridSTA_all{k}, fieldName);

            if all(isnan(vals))
                validField = false;
                break;
            end
        end

        if validField
            musclesToPlot{end+1} = fieldName; 
        end
    end
end

if isempty(musclesToPlot)
    error('No valid muscle fields found to plot.');
end

disp('Muscles that will be plotted:');
disp(musclesToPlot(:));

%% -------------------- CORRECT PAD-SIDE UTAH LAYOUT ----------------------
% Correct electrode layout viewed from pad side.
%
% Top row:
%   9 19 29 39 49 59 69 79
%
% Bottom row:
%   18 28 38 48 58 68 78 88
%
% NaN entries are the four gray non-electrode corner pads.

electrodeLayout = [ ...
    NaN   9   19   29   39   49   59   69   79  NaN;
      1  10   20   30   40   50   60   70   80   89;
      2  11   21   31   41   51   61   71   81   90;
      3  12   22   32   42   52   62   72   82   91;
      4  13   23   33   43   53   63   73   83   92;
      5  14   24   34   44   54   64   74   84   93;
      6  15   25   35   45   55   65   75   85   94;
      7  16   26   36   46   56   66   76   86   95;
      8  17   27   37   47   57   67   77   87   96;
    NaN  18   28   38   48   58   68   78   88  NaN];

% Rotate entire heatmap layout by 180 degrees
electrodeLayout = rot90(electrodeLayout, 2);

%% -------------------- COLORMAP ------------------------------------------

cmap = makeBluePurpleColormap(256);

%% -------------------- PLOT ONE FIGURE PER MUSCLE ------------------------

for m = 1:numel(musclesToPlot)

    muscleField = musclesToPlot{m};
    muscleTitle = prettyMuscleName(muscleField);

    %% -------- Compute muscle-specific colorbar limit across conditions ---

    allValsThisMuscle = [];

    for k = 1:numConditions
        vals = getGridSTAValues(gridSTA_all{k}, muscleField);
        vals = vals(isfinite(vals));
        allValsThisMuscle = [allValsThisMuscle vals]; 
    end

    if isempty(allValsThisMuscle)
        warning('No finite values found for %s. Skipping.', muscleField);
        continue;
    end

    switch colorScaleMode

        case 'trueMaxPerMuscle'
            upperCLim = ceil(max(allValsThisMuscle));

        case 'robust99PerMuscle'
            upperCLim = ceil(prctile(allValsThisMuscle, 99));

        otherwise
            error('Unknown colorScaleMode: %s', colorScaleMode);
    end

    if ~isfinite(upperCLim) || upperCLim <= lowerCLim
        upperCLim = minimumUpperCLim;
    end

    clim = [lowerCLim upperCLim];

    %% -------------------- CREATE FIGURE ---------------------------------

    fig = figure('Color', 'w', 'Position', figurePosition);

    t = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    colormap(fig, cmap);

    axList = gobjects(1, numConditions);

    for k = 1:numConditions

        ax = nexttile(t, k);
        axList(k) = ax;
        hold(ax, 'on');

        heatmapMatrix = buildHeatmapMatrix(gridSTA_all{k}, electrodeLayout, muscleField);

        imagesc(ax, heatmapMatrix, 'AlphaData', ~isnan(heatmapMatrix));

        % Gray background for NaN corner cells
        set(ax, 'Color', [0.65 0.65 0.65]);

        caxis(ax, clim);
        axis(ax, 'image');

        % Row 1 of electrodeLayout is displayed at the top.
        set(ax, 'YDir', 'reverse');

        % Slightly extend the right side so the arrow is visible
        xlim(ax, [0.5 10.8]);
        ylim(ax, [0.5 10.5]);

        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'on';
        ax.LineWidth = 1.5;
        ax.FontName = fontName;

        % Draw grid lines
        drawGridLines(ax, 10, 10);

        % Add electrode numbers
        addElectrodeNumbers(ax, electrodeLayout, numberFontSize, fontName);

        % Add right-pointing arrow
        if showArrow
            q = quiver(ax, 10.18, 5.5, 0.48, 0, 0, ...
                'Color', arrowColor, ...
                'LineWidth', arrowLineWidth, ...
                'MaxHeadSize', 2.5);
            q.Clipping = 'off';
        end

        % Condition title
        title(ax, conditionLabels{k}, ...
            'FontSize', subplotTitleFontSize, ...
            'FontWeight', 'normal', ...
            'FontName', fontName);

        hold(ax, 'off');
    end

    %% -------------------- SHARED COLORBAR -------------------------------

    cb = colorbar(axList(end));
    cb.Layout.Tile = 'east';
    cb.FontName = fontName;
    cb.FontSize = colorbarFontSize;
    cb.LineWidth = 1.5;
    cb.Label.String = 'SD above baseline';
    cb.Label.FontSize = labelFontSize;
    cb.Label.FontName = fontName;
    cb.TicksMode = 'auto';

    %% -------------------- MAIN TITLE ------------------------------------

    title(t, sprintf('%s | %s | StTA response maps', animalLabel, muscleTitle), ...
        'FontSize', mainTitleFontSize, ...
        'FontWeight', 'normal', ...
        'FontName', fontName);

    %% -------------------- SAVE FIGURE -----------------------------------

    if saveFigures

        safeAnimal = regexprep(animalLabel, '[^\w]', '_');
        safeName = regexprep(muscleField, '[^\w]', '_');

        if savePNG
            pngFile = fullfile(saveFolder, sprintf('%s_%s_ThreeConditions.png', ...
                safeAnimal, safeName));
            exportgraphics(fig, pngFile, 'Resolution', 300);
        end

        if savePDF
            pdfFile = fullfile(saveFolder, sprintf('%s_%s_ThreeConditions.pdf', ...
                safeAnimal, safeName));
            exportgraphics(fig, pdfFile, 'ContentType', 'vector');
        end
    end
end

disp('Done plotting three-condition StTA heatmaps.');

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function vals = getGridSTAValues(gridSTA, fieldName)
% Return a 1 x N vector of numeric values from gridSTA.(fieldName).
% Non-numeric, empty, non-scalar, or non-finite entries are returned as NaN.

    vals = nan(1, numel(gridSTA));

    for e = 1:numel(gridSTA)

        if ~isfield(gridSTA(e), fieldName)
            continue;
        end

        value = gridSTA(e).(fieldName);

        if isnumeric(value) && isscalar(value) && isfinite(value)
            vals(e) = value;
        end
    end
end

function heatmapMatrix = buildHeatmapMatrix(gridSTA, electrodeLayout, muscleField)
% Build 10 x 10 heatmap matrix using the physical Utah electrode layout.

    heatmapMatrix = nan(size(electrodeLayout));

    for r = 1:size(electrodeLayout, 1)
        for c = 1:size(electrodeLayout, 2)

            eID = electrodeLayout(r, c);

            if ~isnan(eID)

                if eID > numel(gridSTA)
                    heatmapMatrix(r, c) = NaN;
                    continue;
                end

                if ~isfield(gridSTA(eID), muscleField)
                    heatmapMatrix(r, c) = NaN;
                    continue;
                end

                value = gridSTA(eID).(muscleField);

                if isnumeric(value) && isscalar(value) && isfinite(value)
                    heatmapMatrix(r, c) = value;
                else
                    heatmapMatrix(r, c) = NaN;
                end
            end
        end
    end
end

function addElectrodeNumbers(ax, electrodeLayout, numberFontSize, fontName)
% Add electrode number labels inside each non-corner cell.

    for r = 1:size(electrodeLayout, 1)
        for c = 1:size(electrodeLayout, 2)

            eID = electrodeLayout(r, c);

            if ~isnan(eID)
                text(ax, c, r, sprintf('%d', eID), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', numberFontSize, ...
                    'FontName', fontName, ...
                    'FontWeight', 'bold', ...
                    'Color', [0.15 0.15 0.15]);
            end
        end
    end
end

function cmap = makeBluePurpleColormap(n)
% Light gray/blue to purple colormap similar to the example figure.

    anchorColors = [ ...
        0.95 0.97 0.98;   % very light gray-blue
        0.82 0.89 0.94;   % light blue
        0.64 0.75 0.87;   % medium blue
        0.52 0.58 0.78;   % blue-purple
        0.58 0.34 0.68;   % purple
        0.48 0.05 0.52;   % strong purple
        0.27 0.00 0.32];  % dark purple

    x = linspace(0, 1, size(anchorColors, 1));
    xi = linspace(0, 1, n);

    cmap = interp1(x, anchorColors, xi, 'linear');
end

function drawGridLines(ax, nRows, nCols)
% Draw gray grid lines around all cells.

    gridColor = [0.50 0.50 0.50];
    gridWidth = 1.5;

    for x = 0.5:1:(nCols + 0.5)
        line(ax, [x x], [0.5 nRows + 0.5], ...
            'Color', gridColor, ...
            'LineWidth', gridWidth);
    end

    for y = 0.5:1:(nRows + 0.5)
        line(ax, [0.5 nCols + 0.5], [y y], ...
            'Color', gridColor, ...
            'LineWidth', gridWidth);
    end
end

function prettyName = prettyMuscleName(fieldName)
% Convert field names such as LeftHyoglossus to Left Hyoglossus.

    prettyName = regexprep(fieldName, '([a-z])([A-Z])', '$1 $2');
end