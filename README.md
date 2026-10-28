# PATCH 2.0

MATLAB processing resources for the **Platform for Accurate Tracking of Children’s Health (PATCH)**, a wearable platform for collecting photoplethysmography (PPG) and acceleration signals.


## Repository contents

| File | Purpose |
| --- | --- |
| [HR_spect_smooth_width_V9.m](HR_spect_smooth_width_V9.m) | Main heart-rate processing and export script. |
| [participant_parameters_V6.m](participant_parameters_V6.m) | Participant identifiers, filenames, timing, and analysis settings. |
| [PPG_flatten_V8.m](PPG_flatten_V8.m) | PPG preprocessing and time-domain heart-rate estimation. |
| [accel_overlay_V9.m](accel_overlay_V9.m) | PATCH–Actiheart acceleration alignment and comparison. |
| [activity_lables.m](activity_lables.m) | Activity-log import and interval conversion. The filename spelling is intentional. |
| [accel_V5.m](accel_V5.m) | Earlier acceleration-spectrum and step-estimation workflow. |
| [HR_corr_accel_V5.m](HR_corr_accel_V5.m) | Earlier acceleration-based spectral correction workflow. |
| [PATCH v2.3.7 bill of materials.pdf](PATCH%20v2.3.7%20bill%20of%20materials.pdf) | Hardware bill of materials. |

The repository contains seven MATLAB scripts and one PDF. Participant recordings and activity-log workbooks must be supplied separately.
## Requirements and preparation

Use MATLAB with access to the signal-processing functions used by the scripts, including `pspectrum`, `lowpass`, `findpeaks`, and `xcorr`. Signal Processing Toolbox is therefore needed. Confirm availability of the legacy `nanmean` and `nanstd` functions in your installation. The supported MATLAB release and complete toolbox requirements have not been established by an execution test.

For an initial run:

1. Download or clone this repository.
2. Work in a separate analysis directory containing the scripts and one participant’s input files.
3. Preserve the original recordings and an unchanged copy of the code.
4. Record your MATLAB release, installed toolboxes, repository commit, and any local edits.
5. Confirm the input formats below before processing.

The scripts share variables through the MATLAB workspace. Avoid treating them as independent functions with explicit input and output arguments.

## Input files

### Participant filenames and settings

Edit `participant_parameters_V6.m` for each recording. Its filename conventions are:

| Input | Filename |
| --- | --- |
| PATCH recording | `<ID>_TEST.TXT` |
| Actiheart heart rate | `<ID>_hr.txt` |
| Actiheart acceleration | `<ID>_accel.txt` |

| Setting | Shipped value | Meaning |
| --- | --- | --- |
| `participant_num` | `'1106'` | Identifier used in filenames. |
| `trial_start` | `550` | Initial offset, in seconds. |
| `trial_length` | `2615` | Trial timing parameter; inspect its use before interpreting it as a crop duration. |
| `fs_PPG` | `86.8` | Assumed PATCH sampling frequency in Hz. |
| `HR_start` | `80` | Initial heart rate in bpm; zero requests automatic initialization. |
| `t_Actiheart_shift` | `170` | Initial inter-device offset, in seconds. |
| `time_resolution` | `20` | Spectrogram time resolution. |
| `average_time` | `5` | Acceleration RMS averaging window, in seconds. |

These are participant-specific example settings. The comment suggesting 30-second spectrogram resolution differs from the assigned 20-second value. `specify_start` uses a fixed index conversion; verify its mapping to the actual frequency vector if modifying the spectral configuration.

### PATCH PPG data

`PPG_flatten_V8.m` reads comma-separated data with one header line. It uses numeric column 1 for PPG, starting at sample `1 + round(trial_start * fs_PPG)` and continuing to the file’s end.

Its preprocessing sequence is:

1. Subtract a moving baseline over approximately 0.5 seconds.
2. Apply a 3 Hz IIR low-pass filter with steepness 0.65.
3. Divide by a moving maximum of the absolute baseline-corrected signal.
4. Replace infinite values with missing values and interpolate using splines.
5. Detect peaks above a normalized threshold of 0.5.
6. Calculate time-domain heart rate from successive peak intervals, replace moving-median outliers, and apply a moving mean.

The resulting variables include `PPG_flat`, `HR_flat`, and peak timestamps. The final interval is padded with a one-second value. Window lengths passed without timestamps are counts of observations; for example, smoothing over `time_resolution` peak-derived observations does not guarantee a window of that many seconds.

### Acceleration inputs

`accel_overlay_V9.m` uses PATCH numeric columns 2–4 for the three axes and applies the conversion `(value / 1024) * 25.2 - 16`. It reads tab-separated Actiheart acceleration with six header lines, assumes 50 Hz sampling, and divides vector magnitude by 9.81.

The script computes first differences of vector magnitude, interpolates PATCH differences onto Actiheart timestamps, and estimates an additional delay using cross-correlation over ±300 seconds. It then reconstructs both signals through cumulative sums and removes a four-second moving baseline.

It estimates a scaling factor from moving RMS values and produces acceleration comparison plots. Variables named `ENMO` and `ENMO_actiheart` retain signed values in this implementation. Scaling is applied in comparison plots; the workspace `ENMO` vector itself remains unscaled.

Check both alignment and units before interpreting outputs. The PATCH and Actiheart recordings are cropped using their respective assumed sample rates. The first differences are not divided by sample interval, so they should not automatically be interpreted as physical jerk in g/s.

### Activity log

Supply `<ID>_time log.xlsx` for `activity_lables.m`. The workbook is read into eight variables in this order:

| Position | Imported variable | Content |
| --- | --- | --- |
| 1 | `labels` | Activity names. |
| 2–3 | `acti_start`, `acti_end` | Activity start and end times. |
| 4–5 | `flag_start`, `flag_end` | Flag start and end times. |
| 6 | `flag` | Flag descriptions. |
| 7–8 | `ppg_begin`, `ppg_end` | Additional PPG timing fields. |

The first activity start defines the time origin. Activity intervals are converted to rounded elapsed seconds; flag intervals are converted when their imported type is not a cell array. `ppg_begin` and `ppg_end` are imported but not used by this helper. No example workbook is included, so confirm Excel time types and missing-value handling against a representative file.

## Run the main workflow

After configuring inputs, set MATLAB’s current folder to the analysis directory and run:

```matlab
HR_spect_smooth_width_V9
```

It calls `accel_overlay_V9`, then `PPG_flatten_V8`, and later `activity_lables`.

Spectral heart rate uses 0.7–3.5 Hz and peak tracking.

### Written outputs

| Filename | Content |
| --- | --- |
| `<ID>spect.png` | Spectrogram and heart-rate overlay. |
| `<ID>_metrics.xlsx` | Summary metrics. |
| `<ID>_PATCH_HR.csv` | Times, activities, flags, and paired heart rates. |
| `<ID>_PATCH_ENMO.csv` | Times, activities, flags, and acceleration traces. |
| `<ID>_Self_Consistency.csv` | Time, self-consistency, and width variability. |

Check CSV structure: these exports can contain a single table variable holding a string matrix.

### Interpretation checks

- `Actiheart_RMSE` calculates error standard deviation, not RMSE.
- `Actiheart_MAE_percent` is a fraction; multiply by 100 for percent.
- Accuracy uses an absolute error below 5 bpm; self-consistency uses a difference below 10 bpm.
- Self-consistency windows advance by peak index, despite comments describing one-second steps. Preallocation can leave unused zero rows.
- Reference zeros can become `NaN`; inspect denominators and valid pairs before reporting metrics.

## Earlier analysis scripts

### `accel_V5.m`

This script expects `A` and other variables to exist in the workspace and calls the absent `participant_parameters_V5`. It crops numeric acceleration columns 2–4 between the configured start and an endpoint derived from `trial_length`, computes vector magnitude, subtracts its mean, and generates a 0–10 Hz spectrogram.

Separate spectral thresholds generate step-related and motion-related arrays. Dominant step frequency is integrated into an estimated cumulative step count and plotted. Unlike `accel_overlay_V9.m`, this script does not apply the same conversion from ADC values to g. Treat it as an earlier analysis branch requiring review of dependencies, units, and crop semantics before use.

### `HR_corr_accel_V5.m`

This script calls the absent `HR_spect_smooth_width_V5`, invokes `accel_V5`, divides the PPG spectrum by thresholded acceleration power, and tracks corrected spectral peaks. It also references `HR_polar_interp` and the absent `width_confidence` helper.

The mixed dependencies and reference-device assumptions require reconciliation before this branch can run. Its existence does not mean acceleration-based spectral division is active in the main workflow.

