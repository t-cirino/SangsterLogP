# logP Estimation and Ionisation Classification

## Overview
This script estimates logP, classifies compounds by ionisation state,
and adjusts logP values for ionised compounds using experimental logD data.

### Classification scheme
| Class | Criterion       | Output       |
|-------|-----------------|--------------|
| N     | ΔDP ≤ 0.3       | `N.csv`      |
| MI    | 0.3 < ΔDP ≤ 1.0 | `MI.csv`     |
| HI    | ΔDP > 1.0       | `HI.csv`     |

where **ΔDP = logP − logD_pred_at_exp_pH**  
and **logP = max(predicted logD across pH 0–14)**

For MI and HI compounds, adjusted logP is computed as:

    logP_adj = logD_exp + ΔDP

## Requirements
- Python ≥ 3.7  
- Standard library only (`csv`, `os`, `sys`) — no additional packages needed

## Usage
```bash
python script.py <input_file.csv>
```
Example:
```bash
python script.py example.csv
```
Output files (`N.csv`, `MI.csv`, `HI.csv`) are written to the same
directory as the input file.

## Input file format
CSV with a header row and the following columns:

| Index | Field              | Type  | Notes                        |
|-------|--------------------|-------|------------------------------|
| 0     | SMILES             | str   |                              |
| 1     | ID                 | str   | unique compound identifier   |
| 2     | Experimental pH    | float | must be in range [0, 14]     |
| 3     | Experimental logD  | float | measured at the exp. pH      |
| 4     | Predicted logD 7.4 | float | predicted at pH 7.4          |
| 5–19  | Predicted logD 0–14| float | predicted at pH 0, 1, ..., 14|

Missing or non-numeric values are skipped with a warning.

## Output file format

**N.csv** — neutral-like compounds (no adjustment needed)
| Column | Description               |
|--------|---------------------------|
| ID     | compound identifier       |
| SMILES | molecular structure       |
| logP   | experimental logD (≈ logP)|

**MI.csv / HI.csv** — ionised compounds (adjusted logP)
| Column   | Description                              |
|----------|------------------------------------------|
| ID       | compound identifier                      |
| exp_pH   | experimental pH                          |
| exp_logD | experimentally measured logD             |
| deltaDP  | ΔDP = logP − logD_pred_at_exp_pH         |
| adj_logP | adjusted logP = logD_exp + ΔDP           |
| opt_pH   | pH at which predicted logD is maximum    |

## Citation
If you use this script in your work, please cite the original article:

> Cirino, T., Caron, G., Ermondi, G., Charochkina, L., & Tetko, I. (2026). SangsterLogP - the largest publicly available dataset of logP values. Scientific Data, Volume, Pages.
> https://doi.org/10.xxxx/xxxxxx

The Zenodo repository DOI is:
> https://doi.org/10.5281/zenodo.19387551

## Licence
This script is released under the [MIT License](https://opensource.org/licenses/MIT).
