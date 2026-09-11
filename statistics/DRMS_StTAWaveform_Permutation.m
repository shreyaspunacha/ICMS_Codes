%% ========================================================================
% d_rms MAP DIFFERENCE WITH PULSE-WAVEFORM CONDITION-LABEL PERMUTATION
% BLOCK-PARALLEL VERSION
%
% Statistical procedure:
%   1. Use pulse-triggered rectified waveforms retained during preprocessing
%      by the pulse-level peak selection threshold 
%      (set to zero, so all the pulses).
%   2. Within each electrode and muscle, pool the selected waveforms from
%      the two conditions and permute condition labels while preserving the
%      original number of waveforms assigned to each condition.
%   3. Average the permuted waveforms within each surrogate condition.
%   4. Calculate StTA sdAbove from each surrogate averaged waveform.
%   5. Build surrogate electrode maps for the two conditions.
%   6. Calculate d_rms, global magnitude shift, and spatial correlation r.
%
% The supplied preprocessing script uses a pulse-level peak threshold of
% 0 baseline SD (stdThreshold = 0) to select waveforms. This is done with
% the intension of generalising the code.
%
% Computational approach:
%   - Process one muscle at a time.
%   - Process permutations in blocks to reduce memory use and parallelize
%     the computation.
%   - Preserve the original number of selected waveforms from each condition
%     at every electrode during label permutation.
%   - Use double precision throughout.
%   - Compute surrogate condition B as the complement of surrogate condition A.
%   - Use deterministic Threefry random-number substreams for reproducibility.
%
% Required input variable in each file:
%   DRMS_Input.SelectedWaveforms{electrode, muscle}
%
% Each cell contains:
%   time samples x selected pulse waveforms
%
% Outputs:
%   Results table, diagnostic figures, CSV file, and MAT file.
%% ========================================================================

clear; clc; close all;

%% -------------------- USER SETTINGS -------------------------------------

conditionA_file  = 'DRMS_Input_AllNB_SelectedWaveforms.mat';
conditionB_file  = 'DRMS_Input_SelectiveNB_SelectedWaveforms.mat';

conditionA_label = 'AllNB'; % Combined NB
conditionB_label = 'SelectiveNB';

saveOutputs  = true;
outputPrefix = [conditionA_label '_vs_' conditionB_label '_dRMS_StTAWaveform_Permutation'];

nPerm = 10000;

numWorkers    = 40;
poolProfile   = 'local';
permBlockSize = 250;   % 10000 / 250 = 40 blocks
randomSeed    = 1;

excludeMuscles = {};   % optional muscle names to exclude

alpha = 0.05;

rng(randomSeed);

%% -------------------- LOAD DATA -----------------------------------------

S_A = load(conditionA_file);
S_B = load(conditionB_file);

if ~isfield(S_A, 'DRMS_Input') || ~isfield(S_B, 'DRMS_Input')
    error('Both input files must contain DRMS_Input.');
end

if ~isfield(S_A.DRMS_Input, 'SelectedWaveforms') || ...
        ~isfield(S_B.DRMS_Input, 'SelectedWaveforms')
    error(['Both input files must contain ' ...
        'DRMS_Input.SelectedWaveforms.']);
end

SelectedWaveformsA = S_A.DRMS_Input.SelectedWaveforms;
SelectedWaveformsB = S_B.DRMS_Input.SelectedWaveforms;

muscleNamesAll = string(S_A.DRMS_Input.emg_muscles(:));
numElectrodes  = S_A.DRMS_Input.numElectrodes;
numMusclesAll  = numel(muscleNamesAll);

baselineSamples_10k = S_A.DRMS_Input.baselineSamples_10k;
analysisSamples_10k = S_A.DRMS_Input.analysisSamples_10k;
totalSamples_10k    = baselineSamples_10k + analysisSamples_10k;

if S_B.DRMS_Input.numElectrodes ~= numElectrodes
    error('Condition A and B have different number of electrodes.');
end
if S_B.DRMS_Input.numMuscles ~= numMusclesAll
    error('Condition A and B have different number of muscles.');
end
if S_B.DRMS_Input.baselineSamples_10k ~= baselineSamples_10k || ...
        S_B.DRMS_Input.analysisSamples_10k ~= analysisSamples_10k
    error('Condition A and B use different baseline or analysis window lengths.');
end
if ~isequal(size(SelectedWaveformsA), [numElectrodes, numMusclesAll]) || ...
        ~isequal(size(SelectedWaveformsB), [numElectrodes, numMusclesAll])
    error(['SelectedWaveforms must be a numElectrodes x numMuscles ' ...
        'cell array in both files.']);
end
if any(cellfun(@(x) ~isempty(x) && size(x,1) ~= totalSamples_10k, ...
        SelectedWaveformsA(:))) || ...
        any(cellfun(@(x) ~isempty(x) && size(x,1) ~= totalSamples_10k, ...
        SelectedWaveformsB(:)))
    error(['Every nonempty SelectedWaveforms cell must have ' ...
        'baselineSamples + analysisSamples rows.']);
end
if isfield(S_A.DRMS_Input, 'stdThreshold') && ...
        isfield(S_B.DRMS_Input, 'stdThreshold') && ...
        S_A.DRMS_Input.stdThreshold ~= S_B.DRMS_Input.stdThreshold
    error('Condition A and B were generated using different SD thresholds.');
end
if isfield(S_B.DRMS_Input, 'emg_muscles')
    muscleNames_B = string(S_B.DRMS_Input.emg_muscles(:));
    if numel(muscleNames_B) == numel(muscleNamesAll) && any(muscleNames_B ~= muscleNamesAll)
        warning('Muscle names/order differ between files. Verify before interpreting.');
    end
end

%% -------------------- EXCLUDE UNWANTED MUSCLES --------------------------

keepMuscle = true(numMusclesAll, 1);
for k = 1:numel(excludeMuscles)
    keepMuscle = keepMuscle & ~strcmpi(muscleNamesAll, excludeMuscles{k});
end

muscleIndices = find(keepMuscle);
muscleNames   = muscleNamesAll(keepMuscle);
numMuscles    = numel(muscleNames);

%% -------------------- INITIALIZE OUTPUT ARRAYS --------------------------

d_rms_obs               = nan(numMuscles, 1);
magnitude_shift         = nan(numMuscles, 1);
spatial_r               = nan(numMuscles, 1);
magnitude_component     = nan(numMuscles, 1);
heterogeneous_component = nan(numMuscles, 1);

p_drms           = nan(numMuscles, 1);
p_magnitude      = nan(numMuscles, 1);
p_spatial_change = nan(numMuscles, 1);

MapA = nan(numElectrodes, numMuscles);
MapB = nan(numElectrodes, numMuscles);

d_null_all   = nan(nPerm, numMuscles);
mag_null_all = nan(nPerm, numMuscles);
r_null_all   = nan(nPerm, numMuscles);

%% -------------------- LOAD SELECTED WAVEFORM CELLS ---------------------
% Input cells are electrode x muscle. Transpose the retained muscle subset
% to muscle x electrode for the analysis below.

fprintf('Loading selected waveforms for all muscles...\n');

AllA_waveforms = SelectedWaveformsA(:, muscleIndices).';
AllB_waveforms = SelectedWaveformsB(:, muscleIndices).';

clear SelectedWaveformsA SelectedWaveformsB S_A S_B

%% -------------------- START PARALLEL POOL -------------------------------
% Workers receive only one muscle's selected waveform data at a time.

pool = gcp('nocreate');

if isempty(pool) || pool.NumWorkers ~= numWorkers

    if ~isempty(pool)
        delete(pool);
    end

    clusterObj = parcluster(poolProfile);
    clusterObj.JobStorageLocation = tempdir;
    pool = parpool(clusterObj, numWorkers);
end

fprintf('Parallel pool active with %d workers.\n', pool.NumWorkers);

numBlocks = ceil(nPerm / permBlockSize);

fprintf('Running %d permutations in %d blocks of at most %d permutations.\n', ...
    nPerm, numBlocks, permBlockSize);

%% -------------------- OBSERVED STATS AND BLOCK PERMUTATION --------------

for m = 1:numMuscles

    fprintf('\nProcessing muscle %d/%d: %s\n', ...
        m, numMuscles, muscleNames(m));

    A_wav_m = AllA_waveforms(m, :);
    B_wav_m = AllB_waveforms(m, :);

    %% ---- Observed StTA maps -------------------------------------------

    A_map = nan(numElectrodes, 1);
    B_map = nan(numElectrodes, 1);

    for e = 1:numElectrodes
        A_map(e) = compute_stta_sdabove( ...
            A_wav_m{e}, baselineSamples_10k, analysisSamples_10k);
        B_map(e) = compute_stta_sdabove( ...
            B_wav_m{e}, baselineSamples_10k, analysisSamples_10k);
    end

    MapA(:, m) = A_map;
    MapB(:, m) = B_map;

    %% ---- Observed statistics ------------------------------------------

    d_obs = compute_d_rms(A_map, B_map);
    d_rms_obs(m) = d_obs;

    valid_map = ~isnan(A_map) & ~isnan(B_map);

    if any(valid_map)
        mag_obs = mean(B_map(valid_map)) - mean(A_map(valid_map));
    else
        mag_obs = NaN;
    end

    magnitude_shift(m) = mag_obs;

    r_obs = compute_spatial_r(A_map, B_map);
    spatial_r(m) = r_obs;

    diff_map  = B_map - A_map;
    valid_dec = ~isnan(diff_map);

    if any(valid_dec)
        magnitude_component(m)     = mag_obs^2;
        heterogeneous_component(m) = var(diff_map(valid_dec), 1);
    end

    %% ---- Pool waveforms for permutation for this muscle ---------------

    pooledBase = cell(numElectrodes, 1);
    pooledAnal = cell(numElectrodes, 1);

    totalBaseSum = cell(numElectrodes, 1);
    totalAnalSum = cell(numElectrodes, 1);

    nA_byElectrode = zeros(numElectrodes, 1);
    nB_byElectrode = zeros(numElectrodes, 1);

    for e = 1:numElectrodes

        nA = size(A_wav_m{e}, 2);
        nB = size(B_wav_m{e}, 2);

        nA_byElectrode(e) = nA;
        nB_byElectrode(e) = nB;

        if nA == 0 || nB == 0
            continue;
        end

        pooledWave = [A_wav_m{e}, B_wav_m{e}];

        pooledBase{e} = pooledWave(1:baselineSamples_10k, :);
        pooledAnal{e} = pooledWave( ...
            baselineSamples_10k + 1 : ...
            baselineSamples_10k + analysisSamples_10k, :);

        totalBaseSum{e} = sum(pooledBase{e}, 2);
        totalAnalSum{e} = sum(pooledAnal{e}, 2);
    end

    clear A_wav_m B_wav_m pooledWave

    AllA_waveforms(m, :) = cell(1, numElectrodes);
    AllB_waveforms(m, :) = cell(1, numElectrodes);

    %% ---- Block-parallel permutations ----------------------------------

    d_block   = cell(numBlocks, 1);
    mag_block = cell(numBlocks, 1);
    r_block   = cell(numBlocks, 1);

    parfor blockIdx = 1:numBlocks

        firstPerm = (blockIdx - 1) * permBlockSize + 1;
        lastPerm  = min(blockIdx * permBlockSize, nPerm);
        nThisBlock = lastPerm - firstPerm + 1;

        % Deterministic substream unique to each muscle and permutation block.
        localStream = RandStream('Threefry', 'Seed', randomSeed);
        localStream.Substream = (m - 1) * numBlocks + blockIdx;

        A_perm_map = nan(numElectrodes, nThisBlock);
        B_perm_map = nan(numElectrodes, nThisBlock);

        for e = 1:numElectrodes

            nA = nA_byElectrode(e);
            nB = nB_byElectrode(e);
            nPool = nA + nB;

            if nA == 0 || nB == 0 || nPool < 2
                continue;
            end

            pb = pooledBase{e};
            pa = pooledAnal{e};

            % Generate random condition-label assignments for this electrode and block.
            [~, allPerms] = sort( ...
                rand(localStream, nPool, nThisBlock), 1);

            fakeA_rows = allPerms(1:nA, :);

            % Double-precision assignment/averaging matrix for surrogate A.
            S_A = zeros(nPool, nThisBlock);

            columnOffsets = (0:nThisBlock - 1) * nPool;
            linearIdx = fakeA_rows + columnOffsets;
            S_A(linearIdx) = 1 / nA;

            % Surrogate-A mean baseline and analysis waveforms.
            meanBaseA = pb * S_A;
            meanAnalA = pa * S_A;

            % Surrogate B is the exact complement of surrogate A.
            meanBaseB = (totalBaseSum{e} - nA * meanBaseA) / nB;
            meanAnalB = (totalAnalSum{e} - nA * meanAnalA) / nB;

            % StTA sdAbove for all permutations in this block.
            muA = mean(meanBaseA, 1);
            sdA = std(meanBaseA, 0, 1);
            pkA = max(meanAnalA, [], 1);

            muB = mean(meanBaseB, 1);
            sdB = std(meanBaseB, 0, 1);
            pkB = max(meanAnalB, [], 1);

            sdAboveA = (pkA - muA) ./ sdA;
            sdAboveB = (pkB - muB) ./ sdB;

            sdAboveA(sdA == 0) = NaN;
            sdAboveB(sdB == 0) = NaN;

            A_perm_map(e, :) = sdAboveA;
            B_perm_map(e, :) = sdAboveB;
        end

        % Compute vectorized null statistics using electrodes valid in both maps.
        [d_null_block, mag_null_block, r_null_block] = ...
            compute_null_statistics(A_perm_map, B_perm_map);

        d_block{blockIdx}   = d_null_block(:);
        mag_block{blockIdx} = mag_null_block(:);
        r_block{blockIdx}   = r_null_block(:);
    end

    %% ---- Collect permutation blocks in original order -----------------

    d_null   = vertcat(d_block{:});
    mag_null = vertcat(mag_block{:});
    r_null   = vertcat(r_block{:});

    if numel(d_null) ~= nPerm || ...
            numel(mag_null) ~= nPerm || ...
            numel(r_null) ~= nPerm
        error('Permutation block collection did not produce exactly %d values.', nPerm);
    end

    d_null_all(:, m)   = d_null;
    mag_null_all(:, m) = mag_null;
    r_null_all(:, m)   = r_null;

    %% ---- Observed decomposition sanity check --------------------------

    if ~isnan(d_obs) && ...
            ~isnan(magnitude_component(m)) && ...
            ~isnan(heterogeneous_component(m))

        checkValue = magnitude_component(m) + ...
            heterogeneous_component(m);

        if abs(checkValue - d_obs^2) > 1e-10
            warning('Decomposition mismatch for muscle %d (%s).', ...
                m, muscleNames(m));
        end
    end
end

clear AllA_waveforms AllB_waveforms

%% -------------------- PERMUTATION P-VALUES ------------------------------

for m = 1:numMuscles

    d_obs   = d_rms_obs(m);
    mag_obs = magnitude_shift(m);
    r_obs   = spatial_r(m);

    d_null   = d_null_all(:, m);
    mag_null = mag_null_all(:, m);
    r_null   = r_null_all(:, m);

    % Right-tailed: larger d_rms = larger map difference
    p_drms(m) = (sum(d_null >= d_obs, 'omitnan') + 1) / ...
                (sum(~isnan(d_null)) + 1);

    % Two-tailed: magnitude shift can be positive or negative
    p_magnitude(m) = (sum(abs(mag_null) >= abs(mag_obs), 'omitnan') + 1) / ...
                     (sum(~isnan(mag_null)) + 1);

    % Left-tailed: smaller observed r = reduced spatial similarity
    if ~isnan(r_obs)
        p_spatial_change(m) = (sum(r_null <= r_obs, 'omitnan') + 1) / ...
                              (sum(~isnan(r_null)) + 1);
    end
end

%% -------------------- FDR CORRECTION ------------------------------------

q_drms           = fdr_bh(p_drms);
q_magnitude      = fdr_bh(p_magnitude);
q_spatial_change = fdr_bh(p_spatial_change);

Sig_dRMS_FDR          = q_drms           < alpha;
Sig_Magnitude_FDR     = q_magnitude      < alpha;
Sig_SpatialChange_FDR = q_spatial_change < alpha;

%% -------------------- d_rms² DECOMPOSITION FRACTIONS --------------------

d_rms_sq = d_rms_obs .^ 2;

frac_magnitude     = nan(numMuscles, 1);
frac_heterogeneous = nan(numMuscles, 1);

valid_decomp = d_rms_sq > 0 & ~isnan(d_rms_sq);

frac_magnitude(valid_decomp) = ...
    magnitude_component(valid_decomp) ./ d_rms_sq(valid_decomp);
frac_heterogeneous(valid_decomp) = ...
    heterogeneous_component(valid_decomp) ./ d_rms_sq(valid_decomp);

%% -------------------- INTERPRETATION FRAMEWORK --------------------------

Interpretation = strings(numMuscles, 1);

for m = 1:numMuscles
    sig_d    = Sig_dRMS_FDR(m);
    sig_mag  = Sig_Magnitude_FDR(m);
    sig_spat = Sig_SpatialChange_FDR(m);

    if ~sig_d
        Interpretation(m) = "No significant map difference";
    elseif sig_mag && sig_spat
        Interpretation(m) = "Both magnitude and spatial pattern changed";
    elseif sig_mag && ~sig_spat
        Interpretation(m) = "Global magnitude shift without significant spatial reorganization";
    elseif ~sig_mag && sig_spat
        Interpretation(m) = "Spatial reorganization without global magnitude shift";
    else
        Interpretation(m) = "Distributed change: neither global gain nor spatial r individually significant";
    end
end

%% -------------------- RESULTS TABLE -------------------------------------

Results = table( ...
    muscleNames(:), ...
    d_rms_obs, p_drms, q_drms, Sig_dRMS_FDR, ...
    magnitude_shift, p_magnitude, q_magnitude, Sig_Magnitude_FDR, ...
    spatial_r, p_spatial_change, q_spatial_change, Sig_SpatialChange_FDR, ...
    magnitude_component, heterogeneous_component, ...
    frac_magnitude, frac_heterogeneous, ...
    Interpretation, ...
    'VariableNames', { ...
        'Muscle', ...
        'd_rms', 'p_dRMS', 'q_dRMS_FDR', 'Sig_dRMS_FDR', ...
        'MagnitudeShift_BminusA', 'p_MagnitudeShift', 'q_MagnitudeShift_FDR', 'Sig_MagnitudeShift_FDR', ...
        'SpatialR', 'p_SpatialChange', 'q_SpatialChange_FDR', 'Sig_SpatialChange_FDR', ...
        'MagnitudeComponent', 'HeterogeneousComponent', ...
        'Frac_MagnitudeComponent', 'Frac_HeterogeneousComponent', ...
        'Interpretation'});

disp(Results);

%% -------------------- FIGURE 1: d_rms BAR PLOT --------------------------

fig1 = figure('Color','w','Position',[50 600 900 400]);
bar(d_rms_obs, 'FaceColor', [0.3 0.5 0.8]);
hold on;
add_sig_stars(d_rms_obs, q_drms, alpha, false);
set(gca, 'XTick', 1:numMuscles, 'XTickLabel', muscleNames, ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('d_{rms} of StTA sdAbove maps');
title([conditionA_label ' vs ' conditionB_label ' — Overall Map Difference']);
% Guard against zero or NaN maxima when setting the y-axis limit.
maxD = max(d_rms_obs, [], 'omitnan');
if ~isnan(maxD) && maxD > 0
    ylim([0, maxD * 1.25]);
end
grid on; box on;

%% -------------------- FIGURE 2: MAGNITUDE SHIFT -------------------------

fig2 = figure('Color','w','Position',[50 150 900 400]);
bar(magnitude_shift, 'FaceColor', [0.8 0.4 0.3]);
hold on;
yline(0, 'k--', 'LineWidth', 1.2);
add_sig_stars(magnitude_shift, q_magnitude, alpha, true);
set(gca, 'XTick', 1:numMuscles, 'XTickLabel', muscleNames, ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel(['Mean StTA sdAbove: ' conditionB_label ' - ' conditionA_label]);
title([conditionA_label ' vs ' conditionB_label ' — Global Magnitude Shift']);
grid on; box on;

%% -------------------- FIGURE 3: SPATIAL r -------------------------------

fig3 = figure('Color','w','Position',[1000 600 900 400]);
bar(spatial_r, 'FaceColor', [0.3 0.7 0.4]);
hold on;
yline(0, 'k--', 'LineWidth', 1.2);
add_sig_stars(spatial_r, q_spatial_change, alpha, true);
set(gca, 'XTick', 1:numMuscles, 'XTickLabel', muscleNames, ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Spatial correlation r');
title([conditionA_label ' vs ' conditionB_label ...
    ' — Spatial Pattern Similarity' newline ...
    '(* = significantly reduced spatial similarity, left-tailed permutation test)']);
ylim([-1.1, 1.1]);
grid on; box on;

%% -------------------- FIGURE 4: d_rms² DECOMPOSITION -------------------

fig4 = figure('Color','w','Position',[1000 150 900 400]);
stackData = [frac_magnitude, frac_heterogeneous] * 100;
b4 = bar(stackData, 'stacked');
b4(1).FaceColor = [0.8 0.4 0.3];
b4(2).FaceColor = [0.3 0.7 0.4];
set(gca, 'XTick', 1:numMuscles, 'XTickLabel', muscleNames, ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('% of d_{rms}^2');
title([conditionA_label ' vs ' conditionB_label ' — d_{rms}^2 Decomposition']);
legend({'Magnitude component', 'Heterogeneous component'}, 'Location', 'bestoutside');
ylim([0 110]);
grid on; box on;

%% -------------------- FIGURE 5: INTERPRETATION SUMMARY ------------------

fig5 = figure('Color','w','Position',[550 350 900 450]);
axis off;
title([conditionA_label ' vs ' conditionB_label ' — Interpretation Summary'], ...
    'FontSize', 12, 'FontWeight', 'bold');
subtitle(sprintf('FDR alpha = %.2f | nPerm = %d', alpha, nPerm));

colorMap = containers.Map( ...
    {'No significant map difference', ...
     'Both magnitude and spatial pattern changed', ...
     'Global magnitude shift without significant spatial reorganization', ...
     'Spatial reorganization without global magnitude shift', ...
     'Distributed change: neither global gain nor spatial r individually significant'}, ...
    {[0.75 0.75 0.75], [0.50 0.20 0.70], [0.80 0.40 0.30], ...
     [0.30 0.70 0.40], [0.95 0.75 0.20]});

for m = 1:numMuscles
    label = char(Interpretation(m));
    clr   = colorMap(label);
    yPos  = 1 - m / (numMuscles + 1);
    annotation('textbox', [0.05 yPos 0.90 0.055], ...
        'String', sprintf('%-28s   %s', muscleNames(m), label), ...
        'FitBoxToText', 'off', 'BackgroundColor', clr, 'EdgeColor', 'none', ...
        'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
end

%% -------------------- FIGURE 6: SPATIAL r NULL DISTRIBUTIONS ------------

nCols = ceil(sqrt(numMuscles));
nRows = ceil(numMuscles / nCols);

fig6 = figure('Color','w','Position',[100 100 300*nCols 250*nRows]);

for m = 1:numMuscles
    subplot(nRows, nCols, m);
    r_null_m = r_null_all(:,m);
    r_null_m = r_null_m(~isnan(r_null_m));
    histogram(r_null_m, 50, 'FaceColor', [0.4 0.6 0.8], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.75);
    hold on;
    xline(spatial_r(m), 'r-', 'LineWidth', 2, ...
        'Label', sprintf('obs r=%.2f', spatial_r(m)), ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
    pct5 = prctile(r_null_m, 5);
    xline(pct5, 'k--', 'LineWidth', 1, ...
        'Label', sprintf('5th pct=%.2f', pct5), ...
        'LabelVerticalAlignment', 'top', 'FontSize', 7);
    xlabel('r_{null}', 'FontSize', 8);
    ylabel('Count', 'FontSize', 8);
    if Sig_SpatialChange_FDR(m)
        titleStr = sprintf('%s  [q=%.3f]*', muscleNames(m), q_spatial_change(m));
    else
        titleStr = sprintf('%s  [q=%.3f]', muscleNames(m), q_spatial_change(m));
    end
    title(titleStr, 'FontSize', 8, 'FontWeight', 'bold');
    xlim([-1.1, 1.1]);
    grid on; box on;
end

sgtitle([conditionA_label ' vs ' conditionB_label ...
    ' — Spatial r Null Distributions (left-tailed)' newline ...
    'Red = observed r | Dashed = 5th percentile of null | * = FDR significant'], ...
    'FontSize', 10, 'FontWeight', 'bold');

%% -------------------- SAVE OUTPUTS --------------------------------------

if saveOutputs
    writetable(Results, [outputPrefix '_Results.csv']);
    save([outputPrefix '_Results.mat'], ...
        'Results', 'MapA', 'MapB', 'd_rms_obs', 'p_drms', 'q_drms', ...
        'magnitude_shift', 'p_magnitude', 'q_magnitude', ...
        'spatial_r', 'p_spatial_change', 'q_spatial_change', ...
        'magnitude_component', 'heterogeneous_component', ...
        'frac_magnitude', 'frac_heterogeneous', ...
        'd_null_all', 'mag_null_all', 'r_null_all', ...
        'Interpretation', 'muscleNames', 'conditionA_label', 'conditionB_label', ...
        'numElectrodes', 'nPerm', 'numWorkers', 'permBlockSize', ...
        'randomSeed', 'alpha', '-v7.3');
    saveas(fig1, [outputPrefix '_Figure1_dRMS.png']);
    saveas(fig2, [outputPrefix '_Figure2_MagnitudeShift.png']);
    saveas(fig3, [outputPrefix '_Figure3_SpatialR.png']);
    saveas(fig4, [outputPrefix '_Figure4_Decomposition.png']);
    saveas(fig5, [outputPrefix '_Figure5_InterpretationSummary.png']);
    saveas(fig6, [outputPrefix '_Figure6_SpatialR_NullDiagnostic.png']);
    fprintf('\nSaved outputs with prefix:\n%s\n', outputPrefix);
end

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function sdAbove = compute_stta_sdabove(waveMatrix, baselineSamples, analysisSamples)
    if isempty(waveMatrix)
        sdAbove = NaN;
        return;
    end
    meanRectWave = mean(waveMatrix, 2)';
    sdAbove = compute_stta_sdabove_from_mean(meanRectWave, baselineSamples, analysisSamples);
end

% -------------------------------------------------------------------------

function sdAbove = compute_stta_sdabove_from_mean(meanRectWave, baselineSamples, analysisSamples)
    if isempty(meanRectWave)
        sdAbove = NaN;
        return;
    end
    baselineWave  = meanRectWave(1:baselineSamples);
    muBase = mean(baselineWave);
    sdBase = std(baselineWave);
    analysisWave  = meanRectWave(baselineSamples+1:baselineSamples+analysisSamples);
    peakAmplitude = max(analysisWave);
    if sdBase > 0
        sdAbove = (peakAmplitude - muBase) / sdBase;
    else
        sdAbove = NaN;
    end
end

% -------------------------------------------------------------------------

function [d_null, mag_null, r_null] = ...
    compute_null_statistics(A_perm_map, B_perm_map)
% Vectorized null statistics using only electrodes valid in both maps.
% This matches the paired-valid behavior of compute_d_rms and
% compute_spatial_r.

    validPair = ~isnan(A_perm_map) & ~isnan(B_perm_map);
    nPair = sum(validPair, 1);

    A0 = A_perm_map;
    B0 = B_perm_map;

    A0(~validPair) = 0;
    B0(~validPair) = 0;

    safeCount = max(nPair, 1);

    % d_rms
    diff0 = B0 - A0;
    d_null = sqrt(sum(diff0.^2, 1) ./ safeCount);
    d_null(nPair == 0) = NaN;

    % Magnitude shift
    meanA = sum(A0, 1) ./ safeCount;
    meanB = sum(B0, 1) ./ safeCount;

    mag_null = meanB - meanA;
    mag_null(nPair == 0) = NaN;

    % Spatial correlation
    Ac = A_perm_map - meanA;
    Bc = B_perm_map - meanB;

    Ac(~validPair) = 0;
    Bc(~validPair) = 0;

    numerator = sum(Ac .* Bc, 1);
    denominator = sqrt( ...
        sum(Ac.^2, 1) .* sum(Bc.^2, 1));

    r_null = numerator ./ denominator;
    r_null(nPair <= 2 | denominator == 0) = NaN;
end

% -------------------------------------------------------------------------

function d = compute_d_rms(A_map, B_map)
    A_map = A_map(:);
    B_map = B_map(:);
    valid = ~isnan(A_map) & ~isnan(B_map);
    if sum(valid) == 0
        d = NaN;
    else
        d = sqrt(mean((B_map(valid) - A_map(valid)).^2));
    end
end

% -------------------------------------------------------------------------

function r = compute_spatial_r(A_map, B_map)
    A_map = A_map(:);
    B_map = B_map(:);
    valid = ~isnan(A_map) & ~isnan(B_map);
    if sum(valid) <= 2
        r = NaN;
        return;
    end
    A_valid = A_map(valid);
    B_valid = B_map(valid);
    if std(A_valid) == 0 || std(B_valid) == 0
        r = NaN;
    else
        r = corr(A_valid, B_valid);
    end
end

% -------------------------------------------------------------------------

function q = fdr_bh(p)
    originalSize = size(p);
    p = p(:);
    q = nan(size(p));
    valid   = ~isnan(p);
    p_valid = p(valid);
    if isempty(p_valid)
        q = reshape(q, originalSize);
        return;
    end
    [p_sorted, sortIdx] = sort(p_valid);
    m        = numel(p_sorted);
    q_sorted = p_sorted .* m ./ (1:m)';
    q_sorted = min(q_sorted, 1);
    for k = m-1:-1:1
        q_sorted(k) = min(q_sorted(k), q_sorted(k+1));
    end
    q_valid = nan(size(p_valid));
    q_valid(sortIdx) = q_sorted;
    q(valid) = q_valid;
    q = reshape(q, originalSize);
end

% -------------------------------------------------------------------------

function add_sig_stars(values, q_vals, alpha, signed)
    maxVal = max(abs(values), [], 'omitnan');
    if isempty(maxVal) || isnan(maxVal) || maxVal == 0
        offset = 0.05;
    else
        offset = 0.04 * maxVal;
    end
    for m = 1:numel(values)
        if isnan(values(m)) || isnan(q_vals(m)), continue; end
        if     q_vals(m) < 0.001, starText = '***';
        elseif q_vals(m) < 0.01,  starText = '**';
        elseif q_vals(m) < alpha,  starText = '*';
        else,  continue;
        end
        if signed && values(m) < 0
            yPos = values(m) - offset;
        else
            yPos = values(m) + offset;
        end
        text(m, yPos, starText, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 12);
    end
end
