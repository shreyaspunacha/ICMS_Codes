ICMS Orofacial M1 Analysis

Analysis code associated with a study of tongue-muscle outputs evoked by intracortical microstimulation (ICMS) of orofacial primary motor cortex under different oral sensory conditions. The title of the paper is ``Oral somatosensation shapes lingual muscle maps in macaque orofacial primary motor cortex"

The repository contains code for ICMS delivery, EMG preprocessing, stimulus-triggered averaging (StTA), waveform selection, cortical heatmap generation, and permutation-based comparison of StTA maps across conditions.

## Repository structure

```text
ICMS_Codes/
├── README.md
├── .gitignore
├── stimulation/
│   └── M1_Utah_ICMS_Stimulation.m
├── preprocessing/
│   ├── ICMS_Analysis.m
│   └── ICMS_Waveforms.m
├── figures/
│   └── StTA_ThreeCondition_Heatmaps.m
├── statistics/
│   └── DRMS_StTAWaveform_Permutation.m
└── hpc/
    └── run_DRMS_StTAWaveform_Permutation.sbatch
```

Experimental conditions

The analysis is designed for three conditions:

Control

AllNB (Combined NB)

SelectiveNB

The six tongue-muscle EMG channels included in the shared analysis are:

Right genioglossus

Left genioglossus

Right hyoglossus

Left hyoglossus

Right intrinsic

Left intrinsic

Software requirements

The MATLAB scripts require MATLAB and the relevant toolboxes used by the individual analyses.

Additional dependencies include:

Blackrock/CereStim MATLAB API for ICMS delivery

NPMK for reading Blackrock NSx files

MATLAB Parallel Computing Toolbox for the permutation analysis

Vendor software and hardware-specific libraries are not distributed with this repository.

1. ICMS delivery

File:

stimulation/M1_Utah_ICMS_Stimulation.m

This script configures and delivers ICMS through a CereStim stimulator.

The stimulation protocol implemented in the script includes:

96-electrode Utah array

18 µA stimulation amplitude

Cathodal-first biphasic stimulation

200 µs cathodal phase

65.535 ms interphase interval

200 µs anodal phase

53 µs hardware recovery interval

15 pulses per train

34 trains per electrode

0.5 s inter-electrode interval

Before running, update the path to the locally installed CereStim API.

2. EMG preprocessing and StTA analysis

File:

preprocessing/ICMS_Analysis.m

This script:

Loads the CereStim timing channel.

Identifies stimulation trains.

Maps stimulation order to physical Utah-array electrode IDs.

Computes pulse timing.

Loads EMG recordings.

Band-pass filters continuous EMG from 10–500 Hz.

Full-wave rectifies the filtered EMG.

Extracts EMG windows from -20 ms to +40 ms relative to the anodal-phase onset.

Evaluates pulse-triggered waveforms.

Computes the StTA waveform and peak response magnitude in SD above baseline for each electrode and muscle.

The first 14 pulses of each 15-pulse train are analyzed, giving:

34 trains × 14 pulses = 476 candidate triggers per electrode

Input filenames and local NPMK paths should be updated before use.

3. Pulse-level waveform selection

File:

preprocessing/ICMS_Waveforms.m

This script generates the pulse-level waveform dataset used by the permutation analysis.

For each individual rectified pulse waveform, the pre-trigger baseline mean and SD are calculated and the maximum value in the post-trigger analysis window is evaluated.

The supplied script uses:

stdThreshold = 0, meaning all the waveforms are retained.

In general, a waveform is retained when:

analysisPeak >= baselineMean + stdThreshold × baselineSD

The principal output is:

DRMS_Input.SelectedWaveforms{electrode, muscle}

Each cell contains:

time samples × selected pulse waveforms

Run this script separately for each experimental condition.

4. StTA heatmaps

File:

figures/StTA_ThreeCondition_Heatmaps.m

This script plots cortical StTA response maps for the six tongue muscles across:

Control | All NB | Selective NB

For each muscle, the same color scale is used across all three conditions so that response magnitudes can be compared directly.

The script expects one gridSTA MAT file per condition.

Example generic filenames:

Subject_Control_GridSTA.mat
Subject_AllNB_GridSTA.mat
Subject_SelectiveNB_GridSTA.mat

Each file must contain:

gridSTA

with one entry per physical electrode.

5. Permutation-based map comparison

File:

statistics/DRMS_StTAWaveform_Permutation.m

This script compares StTA maps between two experimental conditions using permutation of selected pulse waveforms.

Within each electrode and muscle:

Selected waveforms from the two conditions are pooled.

Condition labels are randomly reassigned while preserving the original number of selected waveforms in each condition.

Surrogate StTA waveforms are calculated.

Surrogate cortical maps are generated.

Null distributions are computed for:

d_rms

global magnitude shift

spatial correlation r

The default analysis uses:

10,000 permutations
40 parallel workers
250 permutations per block

The tests are:

d_rms: right-tailed permutation test

magnitude shift: two-tailed permutation test

spatial correlation: left-tailed permutation test

Permutation p-values use a +1 correction, and false-discovery-rate correction is applied across muscles.

The script also reports the decomposition:

d_rms² = magnitude component + heterogeneous component

Running the permutation analysis with Slurm

A Slurm submission script is provided:

hpc/run_DRMS_StTAWaveform_Permutation.sbatch

Before running:

Make sure the input MAT filenames in DRMS_StTAWaveform_Permutation.m are correct.

Confirm that MATLAB is available on the cluster.

Edit the module load matlab line if your cluster uses a different MATLAB module name.

Add --account or --partition directives if required by your cluster.

Keep --cpus-per-task consistent with numWorkers in the MATLAB script.

From the repository root:

cd statistics
sbatch ../hpc/run_DRMS_StTAWaveform_Permutation.sbatch

Check job status with:

squeue -u "$USER"

Cancel a job with:

scancel JOB_ID

Typical analysis workflow

ICMS delivery
    ↓
Raw EMG + CereStim timing
    ↓
ICMS_Analysis.m
    ↓
gridSTA / electrode-level StTA responses
    ↓
StTA_ThreeCondition_Heatmaps.m

For waveform-level permutation statistics:

Raw EMG + CereStim timing
    ↓
ICMS_Waveforms.m
    ↓
DRMS_Input.SelectedWaveforms
    ↓
DRMS_StTAWaveform_Permutation.m
    ↓
Permutation statistics, result tables, and diagnostic figures

Data can be obtained from the corresponding author upon request.

The scripts use generic filenames and paths. Users should update local file paths and filenames to point to their own data and software installations.

Reproducibility

The permutation analysis uses deterministic Threefry random-number substreams so that results are reproducible for a fixed random seed and analysis configuration.

Default settings include:

nPerm = 10000;
randomSeed = 1;

