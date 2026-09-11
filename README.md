# ICMS Orofacial M1 Analysis

Analysis code associated with the study:

> **Oral somatosensation shapes lingual muscle maps in macaque orofacial primary motor cortex**

This repository contains code for:

- ICMS delivery
- EMG preprocessing
- Stimulus-triggered averaging (StTA)
- Pulse-level waveform selection
- Cortical StTA heatmap generation
- Permutation-based comparison of StTA maps across conditions

---

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

---

## Experimental conditions

The analysis is designed for three conditions:

- **Control**
- **AllNB (Combined NB)**
- **SelectiveNB**

The six tongue-muscle EMG channels included in the shared analysis are:

- Right genioglossus
- Left genioglossus
- Right hyoglossus
- Left hyoglossus
- Right intrinsic
- Left intrinsic

---

## Software requirements

The MATLAB scripts require MATLAB and the relevant toolboxes used by the individual analyses.

Additional dependencies include:

- **Blackrock/CereStim MATLAB API** for ICMS delivery
- **NPMK** for reading Blackrock NSx files
- **MATLAB Parallel Computing Toolbox** for the permutation analysis

> **Note:** Vendor software and hardware-specific libraries are not distributed with this repository.

---

## 1. ICMS delivery

**File**

```text
stimulation/M1_Utah_ICMS_Stimulation.m
```

This script configures and delivers ICMS through a CereStim stimulator.

### Stimulation protocol

- 96-electrode Utah array
- 18 µA stimulation amplitude
- Cathodal-first biphasic stimulation
- 200 µs cathodal phase
- 65.535 ms interphase interval
- 200 µs anodal phase
- 53 µs hardware recovery interval
- 15 pulses per train
- 34 trains per electrode
- 0.5 s inter-electrode interval

Before running the script, update the path to the locally installed CereStim API.

---

## 2. EMG preprocessing and StTA analysis

**File**

```text
preprocessing/ICMS_Analysis.m
```

This script:

1. Loads the CereStim timing channel.
2. Identifies stimulation trains.
3. Maps stimulation order to physical Utah-array electrode IDs.
4. Computes pulse timing.
5. Loads EMG recordings.
6. Band-pass filters continuous EMG from 10–500 Hz.
7. Full-wave rectifies the filtered EMG.
8. Extracts EMG windows from **−20 ms to +40 ms** relative to anodal-phase onset.
9. Evaluates pulse-triggered waveforms.
10. Computes the StTA waveform and peak response magnitude in SD above baseline for each electrode and muscle.

The first 14 pulses of each 15-pulse train are analyzed:

```text
34 trains × 14 pulses = 476 candidate triggers per electrode
```

Before running the script, update the input filenames and local NPMK path.

---

## 3. Pulse-level waveform selection

**File**

```text
preprocessing/ICMS_Waveforms.m
```

This script generates the pulse-level waveform dataset used by the permutation analysis.

For each individual rectified pulse waveform, the pre-trigger baseline mean and SD are calculated, and the maximum value in the post-trigger analysis window is evaluated.

The supplied script uses:

```matlab
stdThreshold = 0;
```

Therefore, all waveforms are retained.

In general, a waveform is retained when:

```text
analysisPeak >= baselineMean + stdThreshold × baselineSD
```

The principal output is:

```matlab
DRMS_Input.SelectedWaveforms{electrode, muscle}
```

Each cell contains:

```text
time samples × selected pulse waveforms
```

Run this script separately for each experimental condition.

---

## 4. StTA heatmaps

**File**

```text
figures/StTA_ThreeCondition_Heatmaps.m
```

This script plots cortical StTA response maps for the six tongue muscles across:

```text
Control | All NB | Selective NB
```

For each muscle, the same color scale is used across all three conditions so that response magnitudes can be compared directly.

The script expects one `gridSTA` MAT file per condition.

### Example filenames

```text
Subject_Control_GridSTA.mat
Subject_AllNB_GridSTA.mat
Subject_SelectiveNB_GridSTA.mat
```

Each file must contain:

```matlab
gridSTA
```

with one entry per physical electrode.

---

## 5. Permutation-based map comparison

**File**

```text
statistics/DRMS_StTAWaveform_Permutation.m
```

This script compares StTA maps between two experimental conditions using permutation of selected pulse waveforms.

### Procedure

Within each electrode and muscle:

1. Selected waveforms from the two conditions are pooled.
2. Condition labels are randomly reassigned while preserving the original number of selected waveforms in each condition.
3. Surrogate StTA waveforms are calculated.
4. Surrogate cortical maps are generated.
5. Null distributions are computed for:
   - `d_rms`
   - Global magnitude shift
   - Spatial correlation `r`

### Default permutation settings

```text
10,000 permutations
40 parallel workers
250 permutations per block
```

### Statistical tests

- **`d_rms`**: right-tailed permutation test
- **Magnitude shift**: two-tailed permutation test
- **Spatial correlation**: left-tailed permutation test

Permutation p-values use a **+1 correction**, and false-discovery-rate correction is applied across muscles.

The script also reports the decomposition:

```text
d_rms² = magnitude component + heterogeneous component
```

---

## Running the permutation analysis with Slurm

A Slurm submission script is provided:

```text
hpc/run_DRMS_StTAWaveform_Permutation.sbatch
```

Before running:

1. Make sure the input MAT filenames in `DRMS_StTAWaveform_Permutation.m` are correct.
2. Confirm that MATLAB is available on the cluster.
3. Edit the `module load matlab` line if your cluster uses a different MATLAB module name.
4. Add `--account` or `--partition` directives if required by your cluster.
5. Keep `--cpus-per-task` consistent with `numWorkers` in the MATLAB script.

From the repository root:

```bash
cd statistics
sbatch ../hpc/run_DRMS_StTAWaveform_Permutation.sbatch
```

Check job status with:

```bash
squeue -u "$USER"
```

Cancel a job with:

```bash
scancel JOB_ID
```

---

## Typical analysis workflow

### StTA map generation

```text
ICMS delivery
    ↓
Raw EMG + CereStim timing
    ↓
ICMS_Analysis.m
    ↓
gridSTA / electrode-level StTA responses
    ↓
StTA_ThreeCondition_Heatmaps.m
```

### Waveform-level permutation statistics

```text
Raw EMG + CereStim timing
    ↓
ICMS_Waveforms.m
    ↓
DRMS_Input.SelectedWaveforms
    ↓
DRMS_StTAWaveform_Permutation.m
    ↓
Permutation statistics, result tables, and diagnostic figures
```

---

## Data availability

Data can be obtained from the corresponding author upon request.

The scripts use generic filenames and paths. Users should update local file paths and filenames to point to their own data and software installations.

---

## Reproducibility

The permutation analysis uses deterministic Threefry random-number substreams so that results are reproducible for a fixed random seed and analysis configuration.

Default settings include:

```matlab
nPerm = 10000;
randomSeed = 1;
```

