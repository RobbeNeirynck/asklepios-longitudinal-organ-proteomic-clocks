#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pydot
from pytetrad.tools.TetradSearch import TetradSearch
from statsmodels.api import OLS, add_constant
from statsmodels.stats.multitest import multipletests

CSV_PATH = "../Results/agegap_sex_adjusted.csv"
OUT_DIR = "../Results/SEM_Network"

B = 500
SEED = 1
CONSENSUS_FRAC = 0.50
SIG_THRESHOLD = 0.05
BETA_SCALE_MAX = 0.2

LAYOUT_OVERLAP = "scale"
LAYOUT_SEP = "+16"

NEATO = "neato"


def parse_args():
    parser = argparse.ArgumentParser(description="Run SEM FGES network analysis.")
    parser.add_argument("--csv", default=CSV_PATH, help="Path to sex-adjusted age-gap CSV.")
    parser.add_argument(
        "--index-col",
        default=None,
        help="Identifier column. If omitted, visit_id or ParticipantID is detected automatically.",
    )
    parser.add_argument("--n-boot", type=int, default=B, help="Number of bootstrap iterations.")
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed.")
    parser.add_argument("--out-dir", default=OUT_DIR, help="Directory for graph outputs.")
    return parser.parse_args()


def setup_graphviz():
    """Use Graphviz from the active conda environment when available."""
    global NEATO
    graphviz_bin = Path(sys.prefix) / "Library" / "bin"
    neato_exe = graphviz_bin / "neato.exe"

    if neato_exe.exists():
        os.environ["PATH"] = str(graphviz_bin) + os.pathsep + os.environ.get("PATH", "")
        NEATO = str(neato_exe)


def read_agegap_data(csv_path, index_col):
    df = pd.read_csv(csv_path)

    if index_col is None:
        for candidate in ("visit_id", "ParticipantID"):
            if candidate in df.columns:
                index_col = candidate
                break

    if index_col is None or index_col not in df.columns:
        available = ", ".join(df.columns[:10])
        raise ValueError(
            "Could not identify an index column. Use --index-col to specify one. "
            f"Available columns include: {available}"
        )

    df = df.set_index(index_col)
    vars_r1 = sorted(c for c in df.columns if c.endswith("_g_R1_res"))
    vars_r2 = sorted(c for c in df.columns if c.endswith("_g_R2_res"))

    if not vars_r1 or not vars_r2:
        raise ValueError("Input CSV must contain *_g_R1_res and *_g_R2_res columns.")

    organs = sorted({c.split("_")[0] for c in vars_r1})
    return df, vars_r1, vars_r2, organs, index_col


def fges_once(sample, vars_r1, vars_r2):
    ts = TetradSearch(sample)
    ts.use_sem_bic(penalty_discount=0.5)

    for var in vars_r1:
        ts.add_to_tier(1, var)
    for var in vars_r2:
        ts.add_to_tier(2, var)

    ts.set_tier_forbidden_within(1, True)
    ts.set_tier_forbidden_within(2, True)

    # Forbid within-organ edges across time points.
    for var_r1 in vars_r1:
        organ_r1 = var_r1.split("_")[0]
        for var_r2 in vars_r2:
            if organ_r1 == var_r2.split("_")[0]:
                ts.set_forbidden(var_r1, var_r2)

    ts.run_fges()
    return pydot.graph_from_dot_data(ts.get_dot())[0]


def bootstrap_consensus_edges(df, vars_r1, vars_r2, n_boot, seed):
    rng = np.random.default_rng(seed)
    edge_counts = {}

    print(f"Starting FGES bootstrap with {n_boot} samples...")
    for i in range(n_boot):
        if (i + 1) % 100 == 0:
            print(f"  Bootstrap {i + 1}/{n_boot}")

        sample = df.sample(frac=1, replace=True, random_state=int(rng.integers(1e9)))
        graph = fges_once(sample, vars_r1, vars_r2)

        for edge in graph.get_edges():
            source = edge.get_source().strip('"')
            target = edge.get_destination().strip('"')

            if source.split("_")[0] != target.split("_")[0]:
                edge_counts[(source, target)] = edge_counts.get((source, target), 0) + 1

    return {edge for edge, count in edge_counts.items() if (count / n_boot) >= CONSENSUS_FRAC}


def filter_edges_with_ols(df, consensus_edges):
    parents = {}
    for source, target in consensus_edges:
        parents.setdefault(target, []).append(source)

    rows = []
    for target, sources in parents.items():
        fit = OLS(df[target], add_constant(df[sources])).fit()
        for source in sources:
            rows.append((source, target, fit.params[source], fit.pvalues[source]))

    if not rows:
        return pd.DataFrame(columns=["from", "to", "beta", "p", "q"])

    stats = pd.DataFrame(rows, columns=["from", "to", "beta", "p"])
    stats["q"] = multipletests(stats["p"], method="fdr_bh")[1]
    return stats[stats["q"] < SIG_THRESHOLD].copy()


def color_scale_linear(beta):
    midpoint = np.array([180, 180, 180])
    target = np.array([70, 130, 180]) if beta < 0 else np.array([178, 34, 34])
    ratio = min(1.0, abs(beta) / BETA_SCALE_MAX)
    color = midpoint + ratio * (target - midpoint)
    return f"#{int(color[0]):02X}{int(color[1]):02X}{int(color[2]):02X}"


def build_graph(organs, stats_sig, seed):
    graph = pydot.Dot(
        graph_type="digraph",
        layout="neato",
        splines="polyline",
        overlap=LAYOUT_OVERLAP,
        sep=LAYOUT_SEP,
        start=str(seed),
    )

    for organ in organs:
        graph.add_node(
            pydot.Node(
                organ,
                shape="circle",
                style="filled",
                fillcolor="#FFFFFF",
                fontname="Arial",
                fontsize="11",
                penwidth="1.2",
                margin="0.05",
            )
        )

    for _, row in stats_sig.iterrows():
        source = row["from"].split("_")[0]
        target = row["to"].split("_")[0]
        beta_abs = abs(row.beta)

        graph.add_edge(
            pydot.Edge(
                source,
                target,
                color=color_scale_linear(row.beta),
                penwidth=f"{0.8 + 4.0 * (beta_abs / BETA_SCALE_MAX):.2f}",
                arrowsize="0.8",
            )
        )

    return graph


def write_graph(graph, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    dot_file = out_dir / "final_network_fixed.dot"
    png_file = out_dir / "final_network.png"
    svg_file = out_dir / "final_network.svg"

    graph.write_raw(str(dot_file))

    subprocess.run([NEATO, "-Gdpi=600", "-Tpng", str(dot_file), "-o", str(png_file)], check=True)
    subprocess.run([NEATO, "-Tsvg", str(dot_file), "-o", str(svg_file)], check=True)

    print(f"Final graph saved to {png_file} and {svg_file}")


def main():
    args = parse_args()
    setup_graphviz()

    csv_path = Path(args.csv)
    out_dir = Path(args.out_dir)

    if not csv_path.exists():
        sys.exit(f"Error: input file not found: {csv_path}")

    df, vars_r1, vars_r2, organs, index_col = read_agegap_data(csv_path, args.index_col)
    print(f"Loaded {df.shape[0]} rows from {csv_path} using index column '{index_col}'.")

    consensus_edges = bootstrap_consensus_edges(df, vars_r1, vars_r2, args.n_boot, args.seed)
    print(f"Consensus cross-organ edges: {len(consensus_edges)}")

    stats_sig = filter_edges_with_ols(df, consensus_edges)
    print(f"Significant edges after OLS/FDR filtering: {len(stats_sig)}")

    graph = build_graph(organs, stats_sig, args.seed)
    write_graph(graph, out_dir)


if __name__ == "__main__":
    main()
