#!/usr/bin/env python

from __future__ import annotations
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm

from sklearn.linear_model import Lasso
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import argparse

# Default File Paths
DATA_FILE       = "../Data/Metadata.csv"
PROTEIN_FILE    = "../Data/Protein_Dataset_Randomized.csv"
ORGAN_MAP_FILE  = "../Data/Organ_Seq_Id_map.csv"
OUT_DIR         = "../Results"


def parse_args():
    parser = argparse.ArgumentParser(description="Run Organ Proteomic Aging Clocks analysis.")
    parser.add_argument("--data", default=DATA_FILE, help="Path to Metadata CSV (sample annotations).")
    parser.add_argument("--protein", default=PROTEIN_FILE, help="Path to Protein abundance CSV.")
    parser.add_argument("--organ-map", default=ORGAN_MAP_FILE, help="Path to Organ-SeqID mapping CSV.")
    parser.add_argument("--out-dir", default=OUT_DIR, help="Directory to save results.")
    return parser.parse_args()

def group_kfold_once(groups, k=5, seed=42):
    """Generates a single random k-fold partition"""
    rng = np.random.default_rng(seed)
    unique = np.array(sorted(set(groups)))
    rng.shuffle(unique)

    sizes = np.full(k, len(unique) // k)
    sizes[: len(unique) % k] += 1

    cur, folds = 0, []
    for s in sizes:
        folds.append(unique[cur : cur + s])
        cur += s

    for target_groups in folds:
        is_test = np.isin(groups, target_groups)
        yield np.flatnonzero(~is_test), np.flatnonzero(is_test)


def train_bagged_lasso(X, y, groups, n_boot, inner_k, seed, n_jobs):
    """Fits n_boot LASSO models with inner-CV alpha selection."""
    alphas = np.logspace(-5, 1, 50)
    inner_cv = GroupKFold(inner_k)

    def _fit(b):
        rng = np.random.default_rng(seed + b)
        idx = rng.choice(len(X), len(X), replace=True)
        Xb, yb, gb = X[idx], y[idx], groups[idx]

        # Inner CV 
        pipe = Pipeline([("s", StandardScaler()), ("l", Lasso())])
        gs = GridSearchCV(
            pipe, {"l__alpha": alphas}, cv=inner_cv,
            scoring="neg_mean_absolute_error", n_jobs=1
        )
        gs.fit(Xb, yb, groups=gb)

        # Final Fit 
        best_alpha = gs.best_params_["l__alpha"]
        final = Pipeline([("s", StandardScaler()), ("l", Lasso(alpha=best_alpha))])
        final.fit(Xb, yb)
        return final, best_alpha

    if n_jobs == 1:
        return [_fit(b) for b in tqdm(range(n_boot), desc="Bootstrap", leave=False)]
    return Parallel(n_jobs=n_jobs)(delayed(_fit)(b) for b in tqdm(range(n_boot), desc="Bootstrap", leave=False))


def ensemble_predict(model_tuples, X):
    """Averages predictions from a list of models."""
    return np.column_stack([m.predict(X) for m, _ in model_tuples]).mean(axis=1)


def run_organ(organ, seqids, df, groups, out_path):
    """Runs the full nested cross-validation pipeline for a single organ."""
    # Setup output directories and identify features
    od = out_path / organ
    od.mkdir(parents=True, exist_ok=True) # Simplified since plots subdir is gone

    feats = [c for c in seqids if c in df.columns]
    if not feats:
        print(f"Skipping {organ}: No matching proteins found.")
        return

    print(f"--- Processing {organ} ({len(feats)} proteins) ---")
    X = df[feats].to_numpy(float)
    y = df["Age"].to_numpy(float)
    coef_path = od / "coefficients.csv"

    # Outer Cross-Validation Loop
    folds = list(group_kfold_once(groups, k=5, seed=42))
    for fold, (tr_ix, te_ix) in tqdm(enumerate(folds), total=len(folds), desc="Outer CV", leave=False):

        # Train bagged models
        models = train_bagged_lasso(
            X[tr_ix], y[tr_ix], groups[tr_ix],
            n_boot=500, inner_k=5, seed=42 + fold, n_jobs=-1
        )

        # Predict on train and test sets
        y_tr_pred = ensemble_predict(models, X[tr_ix])
        y_te_pred = ensemble_predict(models, X[te_ix])

        # Save Predictions
        for name, ix, preds, subset in [
            ("predictions_train.csv", tr_ix, y_tr_pred, "train"),
            ("predictions_test.csv", te_ix, y_te_pred, "test")
        ]:
            out_df = df.iloc[ix][["SampleID", "ParticipantID", "Round"]].copy()
            out_df["Fold"] = fold
            out_df["Set"] = subset
            out_df["TrueAge"] = y[ix]
            out_df["PredictedAge"] = preds
            out_df.to_csv(od / name, mode="a", index=False, header=not (od / name).exists())

        # Save Coefficients
        recs = []
        for b, (model, alpha) in enumerate(models):
            lasso, scaler = model.named_steps["l"], model.named_steps["s"]
            nz = lasso.coef_ != 0
            
            # Features
            if nz.any():
                coef_orig = lasso.coef_[nz] / scaler.scale_[nz]
                for f, cz, cu in zip(np.array(feats)[nz], lasso.coef_[nz], coef_orig):
                    recs.append({"Fold": fold, "Boot": b, "Feat": f, "Z": cz, "U": cu, "Alpha": alpha})
            
            # Intercept
            int_orig = lasso.intercept_ - (lasso.coef_ * scaler.mean_ / scaler.scale_).sum()
            recs.append({"Fold": fold, "Boot": b, "Feat": "_int_", "Z": lasso.intercept_, "U": int_orig, "Alpha": alpha})
        
        pd.DataFrame(recs).to_csv(coef_path, mode="a", index=False, header=not coef_path.exists())

    # Metadata
    with (od / "metadata.json").open("w") as f:
        json.dump({"n_boot": 500, "n_splits": 5, "seed": 42, "n_feats": len(feats)}, f, indent=2)


def main():
    """Main function to run the analysis."""
    args = parse_args()

    out_dir = Path(args.out_dir).resolve()
    print(f"Reading data from {args.data} and {args.protein}...")

    # Load Data
    meta_df = pd.read_csv(args.data)
    prot_df = pd.read_csv(args.protein)
    
    # Merge datasets on SampleID
    df = pd.merge(meta_df, prot_df, on="SampleID", how="inner")

    # Load Organ Map
    om = pd.read_csv(args.organ_map)

    # Group by ParticipantID to prevent leakage
    groups = df["ParticipantID"].to_numpy()

    for organ in tqdm(om["Organ"].unique(), desc="Organs"):
        seqids = om.loc[om["Organ"] == organ, "SeqID"].tolist()
        run_organ(organ, seqids, df, groups, out_path=out_dir)


if __name__ == "__main__":
    main()