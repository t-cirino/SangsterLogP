"""
logP estimation, ionisation classification, and logP adjustment
---------------------------------------------------------------
Classifies compounds into N / MI / HI based on ΔDP = logP - logD_pred,
and adjusts logP for ionised compounds as logP_adj = logD_exp + ΔDP.

Usage:
    python script.py <input_file.csv>

See README.md for full documentation of input/output formats.
"""

import csv
import os
import sys


# ==============================================================================
# CONSTANTS
# ==============================================================================

PH_MIN = 0          # Minimum valid pH on the prediction grid
PH_MAX = 14         # Maximum valid pH on the prediction grid
PH_SPECIAL = 7.4    # Special pH point included alongside integer grid

# ΔDP thresholds for ionisation classification
THRESHOLD_NEUTRAL  = 0.3    # ΔDP ≤ 0.3  → Neutral-like
THRESHOLD_MODERATE = 1.0    # ΔDP ≤ 1.0  → Moderately ionised (else Highly)

# Expected number of columns in the input CSV
EXPECTED_COLUMNS = 20


# ==============================================================================
# INPUT VALIDATION HELPERS
# ==============================================================================

def safe_float(value, field_name, row_id="unknown"):
    """
    Safely convert a raw CSV string to float.

    Parameters
    ----------
    value : str
        Raw string from a CSV cell.
    field_name : str
        Human-readable field name, used in error messages.
    row_id : str
        Compound ID (or line number) for traceability in warnings.

    Returns
    -------
    float
        The parsed numeric value.

    Raises
    ------
    ValueError
        If the string is empty, a recognised missing-value token,
        or cannot be parsed as a float.
    """
    stripped = value.strip()

    # Recognise common placeholders for missing data
    if stripped in ("", "NA", "N/A", "nan", "NaN", "NULL", "null"):
        raise ValueError(
            f"[ID: {row_id}] Missing value in field '{field_name}': '{value}'"
        )

    try:
        return float(stripped)
    except ValueError:
        raise ValueError(
            f"[ID: {row_id}] Non-numeric value in field '{field_name}': '{value}'"
        )


def validate_pH(exp_pH, row_id="unknown"):
    """
    Check that the experimental pH lies within the valid range [0, 14].

    Parameters
    ----------
    exp_pH : float
        The parsed experimental pH value.
    row_id : str
        Compound ID for traceability in warnings.

    Returns
    -------
    bool
        True if pH is within [0, 14], False otherwise.
    """
    if not (PH_MIN <= exp_pH <= PH_MAX):
        print(
            f"  [WARNING] [ID: {row_id}] Experimental pH {exp_pH} is outside "
            f"the valid range [{PH_MIN}, {PH_MAX}]. Skipping row."
        )
        return False
    return True


# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_csv_data(input_file):
    """
    Load, validate, and parse the input CSV into a list of row dictionaries.

    Each row that passes validation is returned as a dict with clearly named
    fields. Rows that fail any validation check are skipped with a printed
    warning so that a single bad row does not abort the entire run.

    Parameters
    ----------
    input_file : str
        Path to the input CSV file.

    Returns
    -------
    list of dict
        Each dict contains:
            'SMILES'    : str             -- molecular SMILES string
            'ID'        : str             -- compound identifier
            'exp_pH'    : float           -- experimental pH (validated in [0,14])
            'exp_logD'  : float           -- experimentally measured logD
            'logd_74'   : float           -- predicted logD at pH 7.4
            'logd_pH'   : list of float   -- predicted logD at pH 0-14 (15 values)

    Raises
    ------
    FileNotFoundError
        If the specified input file does not exist.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: '{input_file}'")

    parsed_rows = []
    skipped     = 0

    with open(input_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)

        # Always skip the header row
        header = next(reader, None)
        if header is None:
            print("  [WARNING] Input file appears to be empty.")
            return []

        for line_num, parts in enumerate(reader, start=2):

            # --- Column count check ------------------------------------------
            if len(parts) < EXPECTED_COLUMNS:
                print(
                    f"  [WARNING] Line {line_num}: Expected {EXPECTED_COLUMNS} "
                    f"columns, got {len(parts)}. Skipping."
                )
                skipped += 1
                continue

            # --- String fields -----------------------------------------------
            smiles      = parts[0].strip()
            compound_id = parts[1].strip()

            # --- Numeric fields with descriptive error handling --------------
            try:
                exp_pH   = safe_float(parts[2], "experimental pH",    compound_id)
                exp_logD = safe_float(parts[3], "experimental logD",   compound_id)
                logd_74  = safe_float(parts[4], "logD at pH 7.4",      compound_id)
                logd_pH  = [
                    safe_float(parts[5 + i], f"logD at pH {i}", compound_id)
                    for i in range(15)          # pH 0 through 14
                ]
            except ValueError as e:
                print(f"  [WARNING] Line {line_num}: {e}. Skipping row.")
                skipped += 1
                continue

            # --- pH range validation -----------------------------------------
            if not validate_pH(exp_pH, compound_id):
                skipped += 1
                continue

            parsed_rows.append({
                "SMILES":   smiles,
                "ID":       compound_id,
                "exp_pH":   exp_pH,
                "exp_logD": exp_logD,
                "logd_74":  logd_74,
                "logd_pH":  logd_pH,
            })

    print(f"  Loaded {len(parsed_rows)} valid rows from '{input_file}' "
          f"({skipped} row(s) skipped).")
    return parsed_rows


# ==============================================================================
# CORE CALCULATIONS
# ==============================================================================

def compute_logP(logd_74, logd_pH):
    """
    Estimate logP as the maximum predicted logD across all pH points.

    Rationale: the neutral form of a compound dominates where logD is
    highest; therefore the maximum predicted logD across pH 0-14 (and
    the special point 7.4) is used as a proxy for logP.

    Parameters
    ----------
    logd_74 : float
        Predicted logD at pH 7.4.
    logd_pH : list of float
        Predicted logD at pH 0, 1, ..., 14 (15 values, index = pH).

    Returns
    -------
    float
        Estimated logP.
    """
    return max(logd_74, *logd_pH)


def map_pH_to_grid(exp_pH):
    """
    Map an experimental pH to the nearest available prediction grid point.

    The grid consists of integer values 0-14 plus the special point 7.4.
    Any pH within 1e-6 of 7.4 maps to 7.4; all others round to the
    nearest integer.

    Parameters
    ----------
    exp_pH : float
        Experimentally measured pH (pre-validated to be in [0, 14]).

    Returns
    -------
    float or int
        7.4 (float) if exp_pH ≈ 7.4, otherwise the nearest integer in [0, 14].
    """
    if abs(exp_pH - PH_SPECIAL) < 1e-6:
        return PH_SPECIAL
    return int(round(exp_pH))


def get_predicted_logd_at_exp_pH(exp_pH, logd_74, logd_pH):
    """
    Retrieve the predicted logD at the grid point closest to exp_pH.

    Parameters
    ----------
    exp_pH : float
        Experimentally measured pH.
    logd_74 : float
        Predicted logD at pH 7.4.
    logd_pH : list of float
        Predicted logD at pH 0-14.

    Returns
    -------
    float
        Predicted logD at the mapped grid pH.
    """
    mapped = map_pH_to_grid(exp_pH)
    if mapped == PH_SPECIAL:
        return logd_74
    return logd_pH[int(mapped)]     # index = pH integer (0-14)


def compute_delta_dp(logP, logd_pred_at_exp_pH):
    """
    Compute the ionisation indicator ΔDP.

        ΔDP = logP − logD_pred_at_exp_pH

    A larger ΔDP indicates greater ionisation at the experimental pH.

    Parameters
    ----------
    logP : float
        Estimated logP (maximum predicted logD).
    logd_pred_at_exp_pH : float
        Predicted logD at the experimental pH grid point.

    Returns
    -------
    float
        ΔDP value.
    """
    return logP - logd_pred_at_exp_pH


def compute_adjusted_logP(exp_logD, delta_dp):
    """
    Compute the adjusted logP for ionised compounds.

    Formula:
        logP_adj = logD_exp + ΔDP

    Rationale:
        For an ionised compound, the raw logP estimate (max logD) does not
        account for the measured ionisation state. Adding ΔDP to the
        experimental logD anchors the estimate to the real measurement
        while preserving the ionisation offset derived from predictions.

    Parameters
    ----------
    exp_logD : float
        Experimentally measured logD at the compound's pH.
    delta_dp : float
        ΔDP = logP − logD_pred_at_exp_pH.

    Returns
    -------
    float
        Adjusted logP, rounded to 2 decimal places.
    """
    return round(exp_logD + delta_dp, 2)


def find_optimal_pH(logP, logd_74, logd_pH):
    """
    Find the pH at which the predicted logD equals the estimated logP.

    When multiple pH points share the maximum logD, their average is
    returned. This pH represents the most neutral predicted state of
    the compound.

    Parameters
    ----------
    logP : float
        Estimated logP (= maximum predicted logD).
    logd_74 : float
        Predicted logD at pH 7.4.
    logd_pH : list of float
        Predicted logD at pH 0-14.

    Returns
    -------
    float or None
        Optimal pH (rounded to 1 decimal), or None if not determinable.
    """
    optimal_pH_values = []

    if logd_74 == logP:
        optimal_pH_values.append(PH_SPECIAL)       # 7.4

    for i, val in enumerate(logd_pH):
        if val == logP:
            optimal_pH_values.append(float(i))     # integer pH 0-14

    if not optimal_pH_values:
        return None

    return round(sum(optimal_pH_values) / len(optimal_pH_values), 1)


# ==============================================================================
# CLASSIFICATION
# ==============================================================================

def classify_and_enrich(data):
    """
    Classify every compound and attach computed values to its row dict.

    For each compound the function:
        1. Estimates logP (max predicted logD across all pH)
        2. Retrieves the predicted logD at the experimental pH grid point
        3. Computes ΔDP = logP − logD_pred_at_exp_pH
        4. Assigns category: 'N', 'MI', or 'HI'
        5. Attaches all computed fields to the row dict in-place

    Parameters
    ----------
    data : list of dict
        Parsed rows from load_csv_data().

    Returns
    -------
    tuple of three lists (neutral, moderate, high)
        Each list contains enriched row dicts for that category.
    """
    neutral  = []   # ΔDP ≤ 0.3
    moderate = []   # 0.3 < ΔDP ≤ 1.0
    high     = []   # ΔDP > 1.0

    for row in data:
        try:
            logP         = compute_logP(row["logd_74"], row["logd_pH"])
            logd_pred    = get_predicted_logd_at_exp_pH(
                               row["exp_pH"], row["logd_74"], row["logd_pH"]
                           )
            delta_dp     = compute_delta_dp(logP, logd_pred)
            opt_pH       = find_optimal_pH(logP, row["logd_74"], row["logd_pH"])

            # Determine category
            if delta_dp <= THRESHOLD_NEUTRAL:
                category = "N"
            elif delta_dp <= THRESHOLD_MODERATE:
                category = "MI"
            else:
                category = "HI"

            # Attach all derived values to the row for downstream writers
            row["logP"]     = logP
            row["delta_dp"] = delta_dp
            row["opt_pH"]   = opt_pH
            row["category"] = category

            if category == "N":
                neutral.append(row)
            elif category == "MI":
                moderate.append(row)
            else:
                high.append(row)

        except (IndexError, ValueError) as e:
            print(
                f"  [WARNING] Could not classify compound "
                f"'{row.get('ID', 'unknown')}': {e}"
            )

    return neutral, moderate, high


# ==============================================================================
# OUTPUT WRITERS
# ==============================================================================

def write_neutral(neutral, output_dir):
    """
    Write N.csv: neutral-like compounds with ID, SMILES, and logP.

    For neutral compounds (ΔDP ≤ 0.3) no adjustment is required.
    logP reported here is the experimental logD, which for neutral
    compounds is a reliable approximation of logP.

    Output columns:
        ID | SMILES | logP

    Parameters
    ----------
    neutral : list of dict
        Enriched row dicts for neutral-like compounds.
    output_dir : str
        Directory where N.csv will be written.
    """
    out_path = os.path.join(output_dir, "N.csv")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "logP"])

        for row in neutral:
            writer.writerow([
                row["ID"],
                row["SMILES"],
                round(row["exp_logD"], 2),  # exp logD ≈ logP for neutral cpds
            ])

    print(f"  Written: {out_path}  ({len(neutral)} neutral compound(s))")


def write_ionised(records, output_dir, category):
    """
    Write MI.csv or HI.csv: ionised compounds with adjusted logP.

    Output columns:
        ID (as original_ID) | exp_pH | exp_logD |
        deltaDP | adj_logP | opt_pH

    Parameters
    ----------
    records : list of dict
        Enriched row dicts for MI or HI compounds.
    output_dir : str
        Directory where the output CSV will be written.
    category : str
        Either 'MI' or 'HI' — determines the output filename.
    """
    out_path = os.path.join(output_dir, f"{category}.csv")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ID",           # original ID
            "exp_pH",
            "exp_logD",
            "deltaDP",
            "adj_logP",
            "opt_pH",
        ])

        for row in records:
            try:
                adj_logP = compute_adjusted_logP(row["exp_logD"], row["delta_dp"])

                writer.writerow([
                    row['ID'],         # traceability: source ID
                    round(row["exp_pH"],   2),
                    round(row["exp_logD"], 2),
                    round(row["delta_dp"], 3),
                    adj_logP,
                    row["opt_pH"],              # None if indeterminate
                ])

            except KeyError as e:
                print(
                    f"  [WARNING] Missing field {e} for compound "
                    f"'{row.get('ID', 'unknown')}'. Row skipped in {category}.csv."
                )

    print(f"  Written: {out_path}  ({len(records)} {category} compound(s))")


def print_summary(neutral, moderate, high):
    """
    Print a classification summary table to stdout.

    Parameters
    ----------
    neutral : list
    moderate : list
    high : list
    """
    total = len(neutral) + len(moderate) + len(high)
    if total == 0:
        print("  No compounds were classified.")
        return

    print("\n--- Classification Summary ---")
    print(f"  Total  : {total}")
    print(f"  N  (neutral, ΔDP ≤ 0.3)           : "
          f"{len(neutral):4d}  ({len(neutral)/total*100:5.1f}%)")
    print(f"  MI (moderate, 0.3 < ΔDP ≤ 1.0)    : "
          f"{len(moderate):4d}  ({len(moderate)/total*100:5.1f}%)")
    print(f"  HI (high, ΔDP > 1.0)               : "
          f"{len(high):4d}  ({len(high)/total*100:5.1f}%)")


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def main(input_file):
    """
    Execute the full classification and adjustment pipeline.

    Steps
    -----
    1. Load and validate all rows from the input CSV.
    2. Classify every compound into N / MI / HI.
    3. Write per-category output CSVs to the same directory as the input.
    4. Print a summary to stdout.

    Parameters
    ----------
    input_file : str
        Path to the input CSV file supplied on the command line.
    """
    output_dir = os.path.dirname(os.path.abspath(input_file))

    print("=" * 60)
    print(f"Input file : {input_file}")
    print(f"Output dir : {output_dir}")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Step 1: Load data
    # ------------------------------------------------------------------
    print("\n[Step 1] Loading data ...")
    data = load_csv_data(input_file)

    if not data:
        print("  No valid data found. Exiting.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 2: Classify
    # ------------------------------------------------------------------
    print("\n[Step 2] Classifying compounds ...")
    neutral, moderate, high = classify_and_enrich(data)

    # ------------------------------------------------------------------
    # Step 3: Write output files
    # ------------------------------------------------------------------
    print("\n[Step 3] Writing output files ...")
    write_neutral(neutral,  output_dir)
    write_ionised(moderate, output_dir, category="MI")
    write_ionised(high,     output_dir, category="HI")

    # ------------------------------------------------------------------
    # Step 4: Summary
    # ------------------------------------------------------------------
    print_summary(neutral, moderate, high)
    print("\nDone.")


# ==============================================================================
# ENTRY POINT
# ==============================================================================

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(
            "Usage: python script.py <input_file.csv>\n"
            "Example: python script.py example.csv"
        )
        sys.exit(1)

    main(sys.argv[1])