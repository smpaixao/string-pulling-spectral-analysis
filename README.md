# String-pulling kinematic analysis

Custom code and source data for the spectral analysis of string-pulling
kinematics reported in:

> Paixão S. et al. *Peripheral mechanoreceptor plasticity and spinal tactile
> integration converge to control dexterous movements.* Cell Reports (2026).

The analysis quantifies the dominant pulling frequency and the rhythmicity of
forelimb movements during string pulling in control and Zic2 neuron-ablated
mice, from hand trajectories tracked with DeepLabCut. It produces the panels of
Supplementary Figure S3C and the accompanying statistics.

## Contents

```
scripts/   string_pulling_AUGMENTED_v2.3_2026-09-17.py   the analysis
data/      frame_rates.xlsx                              frame rate per animal
           Control/, Zic2/                               one .xlsx per animal
results/   augmented_v2.3_2026-09-17/                    output of the script
```

Each animal file holds four columns, `x_preferred`, `y_preferred`,
`x_non_preferred` and `y_non_preferred`, in pixels, normalised to the length of
the string at the start of each video. Trials are separated by blank rows.
Coordinates are raw DeepLabCut output, with y increasing downward; the script
flips the sign internally.

## Running it

Requires Python 3.9 or later with numpy, scipy, pandas, matplotlib and openpyxl.

```
pip install -r requirements.txt
python scripts/string_pulling_AUGMENTED_v2.3_2026-09-17.py
```

The script resolves paths relative to its own location, so it must stay in
`scripts/`. It takes about a minute and writes to
`results/augmented_v<version>_<date>/`, overwriting that folder on a re-run.
Every output file carries the script version in a `script_version` column or in
`run_info.txt`, so results can always be traced back to the code that produced
them.

Re-running the committed version on the committed data reproduces
`results/augmented_v2.3_2026-09-17/` exactly.

## What the analysis does

Pauses are removed with an adaptive velocity threshold, and each continuous
stretch of pulling lasting at least 1.5 s is analysed as a separate segment.
Segments are linearly detrended, high-pass filtered at 1 Hz to remove slow
postural drift, and normalised to unit variance. Power spectral density is
estimated with Welch's method using a fixed 1.5 s window for every animal
regardless of frame rate, then interpolated onto a common 0.05 Hz grid so that
spectra from animals recorded at different frame rates can be averaged.

The dominant frequency is the peak within 1.5 to 8 Hz. Rhythmicity is the share
of power in that band lying within 0.5 Hz of the peak.

Groups are compared with Welch's t-test or the Mann-Whitney U test according to
a Shapiro-Wilk normality test, with Holm correction across the eight
comparisons. Three animals are excluded for performing fewer than 40 pulling
bouts, leaving 11 per group; the exclusions are listed in the script settings.

All analysis parameters are constants at the top of the script, and the version
history in its header records what changed between versions and why.

## Citing

See `CITATION.cff`. Cite the version DOI of the release matching the published
version of record, not the concept DOI.

## License

MIT, see `LICENSE`.
