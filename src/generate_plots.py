import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
# from pypdf import PdfReader, PdfWriter
import scipy.stats as stats
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages

def read_cmd_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot SV types by repeat category and SV length distribution.")
    parser.add_argument('--tsv', required=True, help='Path to input TSV file')
    parser.add_argument('--out', required=True, help='Output dir')
    return parser.parse_args()

def read_tsv(tsv_file):
    """Read and preprocess the TSV file."""
    df = pd.read_csv(tsv_file, sep='\t')
    df['SVLEN'] = pd.to_numeric(df['SVLEN'], errors='coerce')
    df = df.dropna(subset=['SVLEN'])
    return df

def merge_pdfs(output_pdf, input_pdfs):
    try:
        from pypdf import PdfReader, PdfWriter
    except ImportError:
        print("PyPDF not found. Skipping PDF merge.")
        return
    
    writer = PdfWriter()

    for pdf_file in input_pdfs:
        reader = PdfReader(pdf_file)
        for page in reader.pages:
            writer.add_page(page)

    # Save the merged PDF
    with open(output_pdf, "wb") as f_out:
        writer.write(f_out)

def create_svtype_hist_plot(df, output_pdf, sv_order):
    """Generate vertically stacked subplots for SV types with clean formatting (Matplotlib only)."""
    # Define SV type groups and colors
    legend_labels = {
        'NON_REPETITIVE': ['NON_REPETITIVE'],
        'Tandem repeat': ['HOMO', 'STR', 'VNTR', 'TR'],
        'Mobile element': ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']
    }

    classification_type_order = [ct for group in legend_labels.values() for ct in group]

    colors = {
        'NON_REPETITIVE': '#94d0c5',
        'HOMO': '#beb9d8',
        'STR': '#eb8274',
        'VNTR': "#e68eee",
        'TR': '#87b0d2',
        'LTR': '#f4b567',
        'LINE': '#b4d56f',
        'SINE': '#f6cee0',
        'Retroposon': '#d9d9d8',
        'DNA': '#f9f6b7'
    }

    # Group and order data
    grouped = df.groupby(['SVTYPE', 'CLASSIFICATION']).size().reset_index(name='Count')
    grouped['CLASSIFICATION'] = pd.Categorical(
        grouped['CLASSIFICATION'],
        categories=classification_type_order,
        ordered=True
    )
    grouped['SVTYPE'] = pd.Categorical(grouped['SVTYPE'], categories=sv_order, ordered=True)

    # Filter out SV types with no data
    sv_order_with_data = [sv for sv in sv_order if not grouped[grouped['SVTYPE'] == sv].empty]
    num_sv = len(sv_order_with_data)
    if num_sv == 0:
        print("No data available for any SV types.")
        return

    fig, axes = plt.subplots(nrows=num_sv, ncols=1, figsize=(15, 2.5 * num_sv), sharex=True)
    if num_sv == 1:
        axes = [axes]  # make iterable

    for i, sv in enumerate(sv_order_with_data):
        ax = axes[i]
        sub_df = grouped[grouped['SVTYPE'] == sv]

        # Ensure all classification types are present, fill missing counts with 0
        sub_df = (
            sub_df
            .set_index('CLASSIFICATION')
            .reindex(pd.Index(classification_type_order, name='CLASSIFICATION'))
            .fillna({'Count': 0})
            .reset_index()
        )

        x_positions = range(len(sub_df))
        bar_colors = [colors[cls] for cls in sub_df['CLASSIFICATION']]

        ax.bar(x_positions, sub_df['Count'], color=bar_colors)

        # Set category labels for x-axis (only for last subplot)
        if i == num_sv - 1:
            ax.set_xticks(x_positions)
            ax.set_xticklabels(sub_df['CLASSIFICATION'], rotation=45, ha='right', fontsize=10)
        else:
            ax.set_xticks(x_positions)
            ax.set_xticklabels([])

        # Title inside the plot
        ax.text(
            0.95, 0.9, f"{sv}",
            transform=ax.transAxes,
            fontsize=11,
            fontweight='bold',
            va='top',
            ha='right'
        )

        ax.set_ylabel("")
        ax.grid(False)

    axes[-1].set_xlabel("Classification", fontsize=12, fontweight='bold')
    fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical', fontsize=12, fontweight='bold')

    # Custom legend
    legend_handles = []
    for group_name, types in legend_labels.items():
        legend_handles.append(Patch(facecolor='white', edgecolor='white', label=group_name))
        for t in types:
            legend_handles.append(Patch(facecolor=colors[t], label=f"  {t}"))

    fig.legend(
        handles=legend_handles,
        loc='center left',
        bbox_to_anchor=(0.85, 0.5),
        frameon=False
    )

    # Bold group labels in legend
    for text in fig.legends[0].get_texts():
        if not text.get_text().startswith("  "):  # group name
            text.set_weight('bold')
        else:
            text.set_fontstyle('normal')

    plt.tight_layout(rect=[0.07, 0, 0.85, 1])

    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def create_dist_plot(df, output_pdf):
    df['log2_SV_len'] = np.log2(df['SVLEN'])

    def categorize_sv_type(sv_type):
        if sv_type in ['HOMO', 'STR', 'VNTR', 'TR']:
            return 'Tandem repeat'
        elif sv_type in ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']:
            return 'Mobile element'
        else:
            return 'NON_REPETITIVE'

    df['CLASSIFICATION_CATEGORY'] = df['CLASSIFICATION'].apply(categorize_sv_type)

    ordered_classification_types = ['NON_REPETITIVE', 'HOMO', 'STR', 'VNTR', 'TR', 'LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']

    # Use only two colors
    full_color = '#94d0c5'    # green
    partial_color = '#eb8274' # pink

    num_subplots = len(ordered_classification_types)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(15, 1 * num_subplots), sharex=True)
    fig.suptitle("SV length distribution for each classification type", fontsize=14, fontweight='bold')
    if num_subplots == 1:
        axes = [axes]

    for i, classification_type in enumerate(ordered_classification_types):
        ax = axes[i]
        classification_type_df = df[df['CLASSIFICATION'] == classification_type]

        if not classification_type_df.empty:
            x_min = 4
            x_max = df['log2_SV_len'].max()
            x_grid = np.linspace(x_min, x_max, 1000)

            if classification_type_df['CLASSIFICATION_CATEGORY'].iloc[0] == 'Mobile element':
                complete_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Full']
                fragment_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Partial']

                if len(complete_df) > 1:
                    try:
                        kde = stats.gaussian_kde(complete_df['log2_SV_len'])
                        y = kde(x_grid)
                        y /= y.max()
                        ax.fill_between(x_grid, y, color=full_color, alpha=1)
                        ax.plot(x_grid, y, color='grey', linewidth=0.1)
                    except Exception as e:
                        print(f"Error in KDE (Full) for {classification_type}: {e}")
                if len(fragment_df) > 1:
                    try:
                        kde = stats.gaussian_kde(fragment_df['log2_SV_len'])
                        y = kde(x_grid)
                        y /= y.max()
                        ax.fill_between(x_grid, y, color=partial_color, alpha=0.5)
                        ax.plot(x_grid, y, color='grey', linewidth=0.1)
                    except Exception as e:
                        print(f"Error in KDE (Partial) for {classification_type}: {e}")
            else:
                try:
                    kde = stats.gaussian_kde(classification_type_df['log2_SV_len'])
                    y = kde(x_grid)
                    y /= y.max()
                    ax.fill_between(x_grid, y, color=full_color, alpha=1)
                    ax.plot(x_grid, y, color='grey', linewidth=0.1)
                except Exception as e:
                    print(f"Error in KDE for {classification_type}: {e}")

            # Subplot label
            ax.text(0.01, 0.95, classification_type, transform=ax.transAxes,
                    fontsize=10, fontweight='bold', va='top', ha='left')

            # Formatting
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['0', '1'])
            ax.set_xlim(x_min, x_max)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.tick_params(axis='y', which='major', length=5)
            ax.tick_params(axis='y', which='minor', length=2)
            ax.yaxis.set_minor_locator(plt.FixedLocator([0.5]))
            ax.set_ylabel('')

        custom_decimal_ticks = [50, 100, 330, 1000, 6000, 10000]
        custom_log2_ticks = [np.log2(x) for x in custom_decimal_ticks]
        ax.set_xticks(custom_log2_ticks)

        if i < num_subplots - 1:
            ax.set_xlabel('')
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('log$_2$(SV Length)', fontsize=12, fontweight='bold')
            ax.set_xticklabels([str(x) for x in custom_decimal_ticks])
            ax.xaxis.set_major_locator(plt.FixedLocator(custom_log2_ticks))

        for xtick in ax.get_xticks():
            ax.axvline(x=xtick, linestyle='--', color='gray', linewidth=0.5, zorder=1)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_edgecolor('black')

    fig.text(0.1, 0.5, 'Density', va='center', rotation='vertical', fontsize=12, fontweight='bold')
    # Simple legend
    simple_legend = [
        Patch(facecolor=full_color, edgecolor='black', label='Full'),
        Patch(facecolor=partial_color, edgecolor='black', label='Partial')
    ]
    fig.legend(handles=simple_legend, loc='upper right', bbox_to_anchor=(0.99, 0.99), frameon=False)

    plt.subplots_adjust(right=0.9, hspace=0.9)
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def create_svlength_hist_plot(df, output_pdf):
    df['log2_SV_len'] = np.log2(df['SVLEN'])

    def categorize_sv_type(sv_type):
        if sv_type in ['HOMO', 'STR', 'VNTR', 'TR']:
            return 'Tandem repeat'
        elif sv_type in ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']:
            return 'Mobile element'
        else:
            return 'NON_REPETITIVE'

    df['CLASSIFICATION_CATEGORY'] = df['CLASSIFICATION'].apply(categorize_sv_type)

    ordered_classification_types = ['NON_REPETITIVE', 'HOMO', 'STR', 'VNTR', 'TR', 'LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']

    full_color = '#94d0c5'    # green
    partial_color = '#eb8274' # pink

    num_subplots = len(ordered_classification_types)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(15, 1 * num_subplots), sharex=True)
    fig.suptitle("SV length histogram for each classification type", fontsize=14, fontweight='bold')
    if num_subplots == 1:
        axes = [axes]

    bins = np.linspace(4, df['log2_SV_len'].max(), 50)

    for i, classification_type in enumerate(ordered_classification_types):
        ax = axes[i]
        classification_type_df = df[df['CLASSIFICATION'] == classification_type]

        if not classification_type_df.empty:
            if classification_type_df['CLASSIFICATION_CATEGORY'].iloc[0] == 'Mobile element':
                complete_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Full']
                fragment_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Partial']

                if not complete_df.empty:
                    ax.hist(complete_df['log2_SV_len'], bins=bins, color=full_color, alpha=1, label='Full')
                if not fragment_df.empty:
                    ax.hist(fragment_df['log2_SV_len'], bins=bins, color=partial_color, alpha=0.5, label='Partial')
            else:
                ax.hist(classification_type_df['log2_SV_len'], bins=bins, color=full_color, alpha=1)

            # Subplot label
            ax.text(0.01, 0.95, classification_type, transform=ax.transAxes,
                    fontsize=10, fontweight='bold', va='top', ha='left')

            # Formatting
            ax.set_xlim(4, df['log2_SV_len'].max())
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.tick_params(axis='y', which='major', length=5)
            ax.tick_params(axis='y', which='minor', length=2)

        custom_decimal_ticks = [50, 100, 330, 1000, 6000, 10000]
        custom_log2_ticks = [np.log2(x) for x in custom_decimal_ticks]
        ax.set_xticks(custom_log2_ticks)

        if i < num_subplots - 1:
            ax.set_xlabel('')
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('log$_2$(SV Length)', fontsize=12, fontweight='bold')
            ax.set_xticklabels([str(x) for x in custom_decimal_ticks])
            ax.xaxis.set_major_locator(plt.FixedLocator(custom_log2_ticks))

        for xtick in ax.get_xticks():
            ax.axvline(x=xtick, linestyle='--', color='gray', linewidth=0.5, zorder=1)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_edgecolor('black')

    fig.text(0.08, 0.5, 'Count', va='center', rotation='vertical', fontsize=12, fontweight='bold')
    simple_legend = [
        Patch(facecolor=full_color, edgecolor='black', label='Full'),
        Patch(facecolor=partial_color, edgecolor='black', label='Partial')
    ]
    fig.legend(handles=simple_legend, loc='upper right', bbox_to_anchor=(0.99, 0.99), frameon=False)

    plt.subplots_adjust(right=0.9, hspace=0.9)
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def main():
    # Read command-line arguments
    args = read_cmd_args()
    os.makedirs(args.out, exist_ok=True)

    # Read the TSV file
    df = read_tsv(args.tsv)

    output_pdf_0 = os.path.join(args.out, 'distributions.pdf')
    create_dist_plot(df.copy(), output_pdf_0)

    output_pdf_1 = os.path.join(args.out, 'histograms_svlength.pdf')
    create_svlength_hist_plot(df.copy(), output_pdf_1)
    
    sv_order = ['INS', 'DEL', 'INV', 'DUP', 'BND']
    output_pdf_2 = os.path.join(args.out, 'histograms_svtype.pdf')
    create_svtype_hist_plot(df.copy(), output_pdf_2, sv_order)

    output_pdf = os.path.join(args.out, 'plots.pdf')
    merge_pdfs(output_pdf, [output_pdf_0, output_pdf_1, output_pdf_2])

if __name__ == "__main__":
    main()