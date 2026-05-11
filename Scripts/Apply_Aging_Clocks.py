#!/usr/bin/env python
import argparse
import re
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

PROTEIN_ID_RE = re.compile(r"^\d+-\d+$")

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
MODELS_DIR = ROOT_DIR / "Models"
OUT_FILE = ROOT_DIR / "Results" / "organ_age_predictions.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply the trained organ proteomic aging clocks to raw RFU data."
    )
    parser.add_argument("--data", required=True, help="Path to protein abundance CSV with raw RFU values.")
    parser.add_argument("--id-col", default="SampleID", help="Sample identifier column (default: SampleID).")
    parser.add_argument("--out", default=OUT_FILE, help="Output CSV path (default: Results/organ_age_predictions.csv).")
    return parser.parse_args()


def load_clock_files():
    model_paths = sorted(MODELS_DIR.glob("*.pkl"))
    if not model_paths:
        sys.exit(f"\n[ERROR] No .pkl files found in {MODELS_DIR}")

    clocks = []
    for model_path in model_paths:
        clock_data = joblib.load(model_path)
        required_keys = {"organ", "features", "models"}
        missing_keys = required_keys - set(clock_data)
        if missing_keys:
            sys.exit(f"\n[ERROR] Model file {model_path} is missing keys: {sorted(missing_keys)}")
        clocks.append(
            {
                "path": model_path,
                "organ": clock_data["organ"],
                "features": list(clock_data["features"]),
                "models": clock_data["models"],
            }
        )
    return clocks


def warn_about_column_names(user_df, id_col):
    non_protein_columns = [
        c for c in user_df.columns
        if c != id_col and not PROTEIN_ID_RE.match(str(c))
    ]
    if non_protein_columns:
        print(
            f"\n[WARNING] Found {len(non_protein_columns)} column(s) that do not look like "
            "expected protein IDs ('number-number')."
        )
        print("          These columns will be ignored by the clocks.")
        print("          Examples: " + ", ".join(map(str, non_protein_columns[:10])))


def warn_about_input_scale(user_df, diagnostic_features):
    if not diagnostic_features:
        print(
            "\n[CRITICAL WARNING] No numeric protein-like columns were found for scale checking."
        )
        return

    protein_df = user_df[diagnostic_features].apply(pd.to_numeric, errors="coerce")
    values = protein_df.to_numpy(dtype=float).ravel()
    values = values[np.isfinite(values)]

    if values.size == 0:
        print(
            "\n[CRITICAL WARNING] Protein-like columns were found, but they did not contain "
            "numeric values."
        )
        return

    min_val = float(np.min(values))
    median_val = float(np.median(values))
    p95_val = float(np.percentile(values, 95))
    p99_val = float(np.percentile(values, 99))
    max_val = float(np.max(values))
    negative_fraction = float(np.mean(values < 0))

    print(
        "\nInput protein scale check "
        f"(median={median_val:.2f}, p95={p95_val:.2f}, p99={p99_val:.2f}, max={max_val:.2f})."
    )

    if min_val < 0:
        sys.exit(
            "\n[ERROR] Negative protein values detected "
            f"({negative_fraction:.2%} of checked values). Raw RFU values must be >= 0."
        )
    elif median_val < 30 and p99_val < 50:
        print("\n[WARNING] Protein values look like they may already be log2-transformed.")
    elif p95_val < 50:
        print("\n[WARNING] Protein values are unusually small for raw RFU.")
    else:
        print("Input scale check passed.\n")


def add_missing_required_proteins(user_df, required_features):
    missing_features = sorted(set(required_features) - set(user_df.columns))
    if missing_features:
        print(
            f"\n[WARNING] Missing {len(missing_features)} required protein column(s) "
            "across the loaded clocks."
        )
        print("          These proteins were added and imputed as zeros.")
        print("          Examples: " + ", ".join(missing_features[:20]))
        missing_df = pd.DataFrame(0.0, index=user_df.index, columns=missing_features)
        user_df = pd.concat([user_df, missing_df], axis=1)
    return user_df, missing_features


def prepare_clock_input(user_df, clock_features, organ):
    organ_data = user_df[clock_features].apply(pd.to_numeric, errors="coerce")
    nan_count = int(organ_data.isna().sum().sum())
    organ_data = organ_data.fillna(0.0)

    x_new = organ_data.to_numpy(float)
    if np.nanmin(x_new) < 0:
        sys.exit(f"\n[ERROR] Negative values found while processing the {organ} clock.")

    return np.log2(x_new + 1), nan_count


def main():
    args = parse_args()

    print(f"Loading user data from: {args.data}...")
    try:
        user_df = pd.read_csv(args.data)
    except Exception as e:
        sys.exit(f"\n[ERROR] Could not load data: {e}")

    if args.id_col not in user_df.columns:
        sys.exit(f"\n[ERROR] Identifier column '{args.id_col}' not found in the data.")

    clocks = load_clock_files()
    required_features = sorted({feature for clock in clocks for feature in clock["features"]})

    warn_about_column_names(user_df, args.id_col)

    observed_features = sorted(set(required_features) & set(user_df.columns))
    diagnostic_features = observed_features
    if not diagnostic_features:
        diagnostic_features = [
            c for c in user_df.columns
            if c != args.id_col and PROTEIN_ID_RE.match(str(c))
        ]
    warn_about_input_scale(user_df, diagnostic_features)
    user_df, missing_features = add_missing_required_proteins(user_df, required_features)

    print(f"Found {len(clocks)} organ clocks. Beginning inference...\n")

    imputation_warnings = []
    prediction_frames = []

    for clock in clocks:
        organ = clock["organ"]
        clock_features = clock["features"]
        models = clock["models"]

        x_new, nan_count = prepare_clock_input(user_df, clock_features, organ)
        if nan_count:
            imputation_warnings.append(f"{organ}: Imputed {nan_count} missing/non-numeric values as zeros.")

        predictions_matrix = np.column_stack([m.predict(x_new) for m in models])
        final_predicted_age = predictions_matrix.mean(axis=1)

        prediction_frames.append(
            pd.DataFrame(
                {
                    args.id_col: user_df[args.id_col].values,
                    "Organ": organ,
                    "PredictedAge": final_predicted_age,
                }
            )
        )
        print(f"  [SUCCESS] {organ} ages calculated.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results_df = pd.concat(prediction_frames, ignore_index=True)
    results_df.to_csv(out_path, index=False)
    print(f"\nAll clocks applied. Predictions saved to: {out_path}")

    if missing_features or imputation_warnings:
        print("\n" + "="*50)
        print("DATA PROCESSING REPORT")
        print("="*50)
        if missing_features:
            print(
                f"\n[MISSING PROTEINS] {len(missing_features)} required protein column(s) "
                "were absent and imputed as zeros."
            )
            log_file = out_path.with_name(f"{out_path.stem}_missing_proteins.txt")
            with open(log_file, "w") as f:
                f.write(f"Missing {len(missing_features)} required protein columns:\n")
                f.write(", ".join(missing_features) + "\n")
            print(f"  Detailed list of missing proteins saved to: {log_file}")
        if imputation_warnings:
            print("\n[VALUE IMPUTATION] Missing/non-numeric values inside present protein columns were imputed as zeros:")
            for warning in imputation_warnings:
                print(f"  - {warning}")
        print("="*50 + "\n")


if __name__ == "__main__":
    main()
