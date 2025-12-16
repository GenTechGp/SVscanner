#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib
import matplotlib.cm as cm
import numpy as np
from collections import Counter
from textwrap import wrap

TICKS = [0.0, 0.25, 0.5, 0.75, 1.0]

# Desired fixed plot order
CLASS_ORDER = [
    "Repetitive/Tandem",
    "Repetitive/Tandem/HOMO",
    "Repetitive/Tandem/STR",
    "Repetitive/Tandem/VNTR",
    "Repetitive/Tandem/TR",
    "Repetitive/Mobile",
    "Repetitive/Mobile/SINE",
    "Repetitive/Mobile/LINE",
    "Repetitive/Mobile/LTR",
    "Repetitive/Mobile/DNA",
    "Repetitive/Mobile/Retroposon",
    "Repetitive/Mixed",
    "Repetitive/Mixed/HOMO",
    "Repetitive/Mixed/STR",
    "Repetitive/Mixed/VNTR",
    "Repetitive/Mixed/TR",
    "Repetitive/Mixed/SINE",
    "Repetitive/Mixed/LINE",
    "Repetitive/Mixed/LTR",
    "Repetitive/Mixed/DNA",
    "Repetitive/Mixed/Retroposon",
    "NON_REPETITIVE",
]

# Default traceback types to include in the legend even if absent from data
DEFAULT_TRACEBACK_TYPES = [
    "ABC[MK]",
    "ABD[MK]",
    "ABE[ML]",
    "aBG[M]",
    "AbH[K]",
    "AbH[L]",
    "abI",
    "AH[K]",
    "AH[L]",
    "aI",
    "BG[M]",
    "bI",
    "I",
]

def read_cmd_args():
    def positive_float(value):
        ivalue = float(value)
        if ivalue <= 0 or ivalue > 1:
            raise argparse.ArgumentTypeError(f"{value} is an invalid. It must be a value between 0 and 1")
        return ivalue
    
    parser = argparse.ArgumentParser(
        description=(
            "Generate traceback scatter plots from a single TSV containing coverage, "
            "classification, transposition, and traceback."
        )
    )
    parser.add_argument("--traceback", required=True, help="Path to traceback TSV. Lines starting with '#' are treated as comments and parsed as header annotations. "
             "Each data row must be: vcf_id, TRF_TOTAL_SV_COVERAGE, RM_TOTAL_SV_COVERAGE, FINAL_CLASSIFICATION, transposition, Traceback",
    )
    parser.add_argument("--output", required=True, help="Output pdf file path")
    parser.add_argument('--sv_coverage_cutoff', required=False, type=positive_float, default=0.5, help="The soft minimum total sv coverage by repeat elements to be considered repetitive")
    parser.add_argument('--min_total_sv_coverage', required=False, type=positive_float, default=0.75, help="The tight minimum total sv coverage by repeat elements to be considered repetitive")

    return parser.parse_args()

def read_header_comments(path):
    raw_lines = []
    parsed_entries = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                s = line.rstrip('\n')
                raw_lines.append(s)
                content = s[1:].strip()
                if not content:
                    parsed_entries.append(None)
                    continue
                if '\t' in content:
                    parts = content.split('\t', 1)
                elif ' - ' in content:
                    parts = content.split(' - ', 1)
                else:
                    parts = content.split(None, 1)
                if len(parts) == 2:
                    symbol, desc = parts[0].strip(), parts[1].strip()
                    parsed_entries.append((symbol, desc))
                else:
                    parsed_entries.append(None)
            else:
                break
    return raw_lines, parsed_entries

def load_traceback_tsv(path):
    # Read TSV data while skipping explanatory header lines starting with '#'
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
        names=[
            "vcf_id",
            "TRF_TOTAL_SV_COVERAGE",
            "RM_TOTAL_SV_COVERAGE",
            "FINAL_CLASSIFICATION",
            "transposition",
            "Traceback",
        ],
        engine="python",
    )

    # Treat '.' as zero for coverage columns, then convert to numeric
    for col in ["TRF_TOTAL_SV_COVERAGE", "RM_TOTAL_SV_COVERAGE"]:
        df[col] = df[col].replace({'.': '0', None: '0'})
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    # Normalize transposition values (fix ambiguous Series truth-value error)
    # Use fillna + astype + str.strip instead of "or" on a Series
    df["transposition"] = df["transposition"].fillna("").astype(str).str.strip()

    return df

def style_axes(ax, s_cutoff, t_cutoff):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks(TICKS)
    ax.set_yticks(TICKS)
    # Draw cutoff lines from user args
    ax.axhline(s_cutoff, ls="--", color="gray")
    ax.axhline(t_cutoff, ls="--", color="gray")
    ax.axvline(s_cutoff, ls="--", color="gray")
    ax.axvline(t_cutoff, ls="--", color="gray")
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")

def build_distinct_colors(labels):
    labels_sorted = sorted(labels)
    n = len(labels_sorted)

    tab20  = cm.get_cmap('tab20')
    tab20b = cm.get_cmap('tab20b')
    tab20c = cm.get_cmap('tab20c')

    colors = []
    for cmap in [tab20, tab20b, tab20c]:
        for i in range(cmap.N):
            colors.append(cmap(i))
            if len(colors) >= n:
                break
        if len(colors) >= n:
            break

    if len(colors) < n:
        hues = np.linspace(0, 1, n - len(colors), endpoint=False)
        colors.extend([cm.hsv(h) for h in hues])

    return {lab: colors[i] for i, lab in enumerate(labels_sorted)}


def add_legend_outside(ax, traceback_colors, counts, title="Traceback"):
    legend_elements = []
    for lab in traceback_colors.keys():  # preserve insertion order
        count = counts.get(lab, 0)
        label_text = f"{lab} - {count}" if count > 0 else lab
        legend_elements.append(
            Line2D([0], [0], marker="o", color="w", label=label_text,
                   markerfacecolor=traceback_colors[lab], markersize=8)
        )
    ax.legend(
        handles=legend_elements,
        title=title,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=False,
    )

def scatter_points(ax, sub, traceback_colors):
    """
    Scatter points with:
    - color by Traceback
    - marker shape by transposition: triangle (marker='^') for 'Full', circle ('o') otherwise
    """
    if sub.empty:
        return
    # Split by transposition
    full = sub[sub["transposition"].str.lower() == "full"]
    other = sub[sub["transposition"].str.lower() != "full"]

    # Plot 'other' as circles
    if not other.empty:
        colors_other = [traceback_colors.get(t, "#888888") for t in other["Traceback"]]
        ax.scatter(
            other["TRF_TOTAL_SV_COVERAGE"].astype(float),
            other["RM_TOTAL_SV_COVERAGE"].astype(float),
            c=colors_other,
            alpha=0.8,
            s=22,
            edgecolors="none",
            marker='o',
        )

    # Plot 'Full' as small triangles
    if not full.empty:
        colors_full = [traceback_colors.get(t, "#888888") for t in full["Traceback"]]
        ax.scatter(
            full["TRF_TOTAL_SV_COVERAGE"].astype(float),
            full["RM_TOTAL_SV_COVERAGE"].astype(float),
            c=colors_full,
            alpha=0.9,
            s=28,
            edgecolors="none",
            marker='^',
        )

def plot_class_page(df, class_name, traceback_colors, s_cutoff, t_cutoff):
    sub = df[df["FINAL_CLASSIFICATION"] == class_name]
    fig, ax = plt.subplots(figsize=(8.0, 5))  # wider for outside legend

    # Scatter with marker by transposition
    scatter_points(ax, sub, traceback_colors)

    style_axes(ax, s_cutoff, t_cutoff)
    count = len(sub)
    ax.set_title(f"{class_name} ({count})")
    ax.set_xlabel("TRF total SV coverage")
    ax.set_ylabel("RM total SV coverage")

    # Counts per traceback type for this page
    counts = Counter(sub["Traceback"]) if not sub.empty else {}
    add_legend_outside(ax, traceback_colors, counts)
    fig.tight_layout()
    return fig

def plot_header_page(raw_lines, parsed_entries):
    fig, ax = plt.subplots(figsize=(8.5, 11))  # portrait
    ax.axis('off')

    # Title
    ax.text(0.02, 0.98, "Traceback annotations key", fontsize=16, fontweight='bold', va='top', ha='left')

    # Print lines in file order. Prefer parsed_entries when available, otherwise raw without '#'.
    y = 0.94
    for i, raw in enumerate(raw_lines):
        entry = parsed_entries[i]
        if entry is not None:
            sym, desc = entry
            line = f"{sym} - {desc}"
        else:
            line = raw.lstrip('#').strip()
        wrapped = wrap(line, width=95)
        for w in wrapped:
            ax.text(0.02, y, w, fontsize=11, va='top', ha='left')
            y -= 0.03
            if y < 0.05:
                break
        if y < 0.05:
            break

    # Add a note about markers
    ax.text(0.02, max(y-0.02, 0.05),
            "Marker shapes: triangle = transposition 'Full'; circle = other/empty.",
            fontsize=11, va='top', ha='left')

    return fig

def main():
    args = read_cmd_args()

    # Read header comments and mapping for the first page
    raw_lines, parsed_entries = read_header_comments(args.traceback)

    df = load_traceback_tsv(args.traceback)

    # Two-pass: union of default traceback types and types found in data.
    # Preserve insertion order: first defaults in given order, then found types
    found_types = [t for t in df["Traceback"].dropna().unique().tolist()]
    all_types_ordered = DEFAULT_TRACEBACK_TYPES + [t for t in found_types if t not in DEFAULT_TRACEBACK_TYPES]

    # Print found types
    if not all_types_ordered:
        print("No traceback types found in the input.")
    else:
        print("Found traceback types (including defaults, in order):")
        for t in all_types_ordered:
            print(f"- {t}")

    # Build distinct colors in the same order (no sorting)
    traceback_colors_map = build_distinct_colors(all_types_ordered)
    traceback_colors = {lab: traceback_colors_map[lab] for lab in all_types_ordered}  # preserve order

    # User-specified cutoffs
    S = float(args.sv_coverage_cutoff)
    T = float(args.min_total_sv_coverage)

    out_pdf = args.output
    with PdfPages(out_pdf) as pdf:
        # First page: header annotations (original order), with marker note
        header_fig = plot_header_page(raw_lines, parsed_entries)
        pdf.savefig(header_fig, bbox_inches="tight")
        plt.close(header_fig)

        # Subsequent pages: fixed order of classes
        for cls in CLASS_ORDER:
            fig = plot_class_page(df, cls, traceback_colors, S, T)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"Wrote PDF: {out_pdf}")

if __name__ == "__main__":
    # if pdfpages is not found return
    if 'PdfPages' not in globals():
        print("Error: matplotlib PdfPages module not found. Please ensure matplotlib is installed.")
        exit(0)
    main()