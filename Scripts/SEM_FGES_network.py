#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import pydot
import argparse
from statsmodels.api import OLS, add_constant
from statsmodels.stats.multitest import multipletests
from pytetrad.tools.TetradSearch import TetradSearch

# --- Configuration ---
# File and Data Settings
CSV_PATH = "../Results/agegap_sex_adjusted.csv"
INDEX_COL = "ParticipantID"

# Bootstrap & Consensus Settings
B = 500
SEED = 1
CONSENSUS_FRAC = 0.50

# Statistical & Visualization Settings
SIG_THRESHOLD = 0.05

def parse_args():
    parser = argparse.ArgumentParser(description="Run SEM FGES Network Analysis.")
    parser.add_argument("--csv", default=CSV_PATH, help="Path to input CSV file (age gaps).")
    parser.add_argument("--n-boot", type=int, default=B, help="Number of bootstrap iterations.")
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed.")
    return parser.parse_args()

BETA_SCALE_MAX = 0.2

# Graphviz Layout Parameters (for reproducibility)
LAYOUT_SEED = "64652"
LAYOUT_SEP = "+16"
LAYOUT_OVERLAP = "scale"
LAYOUT_K = "0.3"

# Global variable for NEATO executable
NEATO = "neato"

def setup_graphviz():
    """Configures the path to Graphviz binaries."""
    global NEATO
    try:
        graphviz_bin = os.path.join(sys.prefix, "Library", "bin")
        possible_neato = os.path.join(graphviz_bin, "neato.exe")
        if os.path.exists(graphviz_bin) and os.path.exists(possible_neato):
            NEATO = possible_neato
            os.environ["PATH"] = graphviz_bin + os.pathsep + os.environ.get("PATH", "")
    except Exception:
        pass # Fallback to default "neato"

def fges_once(sample: pd.DataFrame, vars_R1: list, vars_R2: list) -> pydot.Graph:
    """Runs a single instance of the FGES algorithm with temporal constraints."""
    ts = TetradSearch(sample)
    ts.use_sem_bic(penalty_discount=0.5)
    
    # Define temporal tiers: R1 (baseline) -> R2 (follow-up)
    for v in vars_R1: ts.add_to_tier(1, v)
    for v in vars_R2: ts.add_to_tier(2, v)
    ts.set_tier_forbidden_within(1, True)
    ts.set_tier_forbidden_within(2, True)
    
    # Forbid self-loops from R1 to R2 for the same organ
    for v1 in vars_R1:
        v2 = v1.replace("_g_R1_res", "_g_R2_res")
        if v1 in sample.columns and v2 in sample.columns:
            ts.set_forbidden(v1, v2)
            
    ts.run_fges()
    return pydot.graph_from_dot_data(ts.get_dot())[0]

def color_scale_linear(beta):
    """Maps a beta value to a color (blue for negative, red for positive)."""
    mid = np.array([180, 180, 180]) # Grey
    ratio = min(1.0, abs(beta) / BETA_SCALE_MAX)
    if beta < 0:
        target = np.array([70, 130, 180]) # Blue
    else:
        target = np.array([178, 34, 34]) # Red
    res = mid + ratio * (target - mid)
    return f"#{int(res[0]):02X}{int(res[1]):02X}{int(res[2]):02X}"

def main():
    args = parse_args()
    csv_path = args.csv
    n_boot = args.n_boot
    seed = args.seed

    setup_graphviz()
    
    # 1. Load and Prepare Data
    if not os.path.exists(csv_path):
        print(f"Error: Input file '{csv_path}' not found.")
        print("Please run the SEM_prepare_agegaps.R script first to generate this file.")
        return

    df = pd.read_csv(csv_path).set_index(INDEX_COL)
    vars_R1 = sorted([c for c in df.columns if c.endswith("_g_R1_res")])
    vars_R2 = sorted([c for c in df.columns if c.endswith("_g_R2_res")])
    organs = sorted({c.split("_")[0] for c in vars_R1})

    print(f"Loaded {df.shape[0]} participants from {csv_path}.")

    # 2. Define and Run Bootstrapped FGES
    print(f"Starting FGES bootstrap with {n_boot} samples...")
    rng = np.random.default_rng(seed)
    edge_counts = {}

    for i in range(n_boot):
        if (i + 1) % 100 == 0:
            print(f"  Bootstrap iteration {i+1}/{n_boot}")
        sample = df.sample(frac=1, replace=True, random_state=int(rng.integers(1e9)))
        g_tmp = fges_once(sample, vars_R1, vars_R2)
        for edge in g_tmp.get_edges():
            u = edge.get_source().strip('"')
            v = edge.get_destination().strip('"')
            edge_counts[(u, v)] = edge_counts.get((u, v), 0) + 1

    consensus_edges = {e for e, c in edge_counts.items() if (c / n_boot) >= CONSENSUS_FRAC}

    # 3. Filter Consensus Edges with OLS Regression
    print("Filtering consensus edges with OLS regression...")
    parents = {}
    for u, v in consensus_edges:
        parents.setdefault(v, []).append(u)

    ols_rows = []
    for to_var, from_vars in parents.items():
        fit = OLS(df[to_var], add_constant(df[from_vars])).fit()
        for frm in from_vars:
            ols_rows.append((frm, to_var, fit.params[frm], fit.pvalues[frm]))

    if not ols_rows:
        print("No significant edges found after OLS filtering.")
        return

    stats = pd.DataFrame(ols_rows, columns=["from", "to", "beta", "p"])
    stats["q"] = multipletests(stats["p"], method="fdr_bh")[1]
    stats_sig = stats[stats["q"] < SIG_THRESHOLD].copy()

    # 5. Generate Graph Layout (Step 4 is the helper function)
    print("Generating graph layout with neato...")
    g_layout = pydot.Dot(
        graph_type="digraph", layout="neato", splines="curved",
        overlap=LAYOUT_OVERLAP, sep=LAYOUT_SEP, start=LAYOUT_SEED, K=LAYOUT_K
    )
    for org in organs:
        g_layout.add_node(pydot.Node(org))
    for _, r in stats_sig.iterrows():
        src = r["from"].split("_")[0]
        dst = r["to"].split("_")[0]
        g_layout.add_edge(pydot.Edge(src, dst))

    layout_dot = g_layout.create_dot(prog=NEATO).decode("utf-8")
    layout_graph = pydot.graph_from_dot_data(layout_dot)[0]

    positions = {}
    for node in layout_graph.get_nodes():
        name = node.get_name().strip('"')
        pos = node.get_pos()
        if pos:
            positions[name] = pos.strip('"')

    # 6. Create and Render the Final Graph
    print("Rendering final graph...")
    g_final = pydot.Dot(graph_type="digraph", layout="neato", splines="polyline")

    for org in organs:
        g_final.add_node(pydot.Node(
            org, shape="circle", style="filled", fillcolor="#FFFFFF",
            fontname="Arial", fontsize="11", penwidth="1.2", margin="0.05",
            pos=positions.get(org), pin="true"
        ))

    for _, r in stats_sig.iterrows():
        src = r["from"].split("_")[0]
        dst = r["to"].split("_")[0]
        beta_abs = abs(r.beta)
        g_final.add_edge(pydot.Edge(
            src, dst,
            color=color_scale_linear(r.beta),
            penwidth=f"{0.8 + 4.0 * (beta_abs / BETA_SCALE_MAX):.2f}",
            arrowsize="0.8"
        ))

    # 7. Write Output Files
    dot_file = "final_network_fixed.dot"
    g_final.write_raw(dot_file)

    output_base = f"final_network_seed{LAYOUT_SEED}"
    png_file = f"{output_base}.png"
    svg_file = f"{output_base}.svg"

    try:
        subprocess.run([NEATO, "-n2", "-Gdpi=600", "-Tpng", dot_file, "-o", png_file], check=True)
        subprocess.run([NEATO, "-n2", "-Tsvg", dot_file, "-o", svg_file], check=True)
        print(f"✅ Final graph saved to {png_file} and {svg_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error generating graph images: {e}")
        print("Ensure 'neato' is installed and in your PATH (part of Graphviz).")

if __name__ == "__main__":
    main()
