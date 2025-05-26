import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pypdf import PdfReader, PdfWriter
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

def create_hist_plot(df, output_pdf, sv_order):
    """Generate vertically stacked subplots for SV types with clean formatting."""
    # Define SV type groups and colors
    legend_labels = {
        'NON_REPETITIVE': ['NON_REPETITIVE'],
        'Tandem repeat': ['HOMO', 'STR', 'TR'],
        'Mobile element': ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']
    }

    classification_type_order = [ct for group in legend_labels.values() for ct in group]

    colors = {
        'NON_REPETITIVE': '#94d0c5',
        'HOMO': '#beb9d8',
        'STR': '#eb8274',
        'TR': '#87b0d2',
        'LTR': '#f4b567',
        'LINE': '#b4d56f',
        'SINE': '#f6cee0',
        'Retroposon': '#d9d9d8',
        'DNA': '#f9f6b7'
    }

    grouped = df.groupby(['SVTYPE', 'CLASSIFICATION']).size().reset_index(name='Count')
    grouped['CLASSIFICATION'] = pd.Categorical(grouped['CLASSIFICATION'], categories=classification_type_order, ordered=True)
    grouped['SVTYPE'] = pd.Categorical(grouped['SVTYPE'], categories=sv_order, ordered=True)

    # Filter out SV types with no data
    sv_order_with_data = [sv for sv in sv_order if not grouped[grouped['SVTYPE'] == sv].empty]

    num_sv = len(sv_order_with_data)
    if num_sv == 0:
        print("No data available for any SV types.")
        return

    sns.set(style="white")  # remove background grid
    fig, axes = plt.subplots(nrows=num_sv, ncols=1, figsize=(15, 2.5 * num_sv), sharex=True)

    if num_sv == 1:
        axes = [axes]  # make iterable

    for i, sv in enumerate(sv_order_with_data):
        ax = axes[i]
        sub_df = grouped[grouped['SVTYPE'] == sv]

        sns.barplot(
            data=sub_df, x='CLASSIFICATION', y='Count', hue='CLASSIFICATION',
            hue_order=classification_type_order, palette=colors,
            dodge=False, ax=ax
        )

        # Move title inside the plot
        ax.text(
            0.95, 0.9, f"{sv}",
            transform=ax.transAxes,
            fontsize=11,
            fontweight='bold',
            va='top',
            ha='right'
        )
        
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(False)  # remove gridlines

        # Remove individual legends
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

    axes[-1].tick_params(axis='x', which='both', bottom=True, labelbottom=True, labelrotation=45, labelsize=10)
    axes[-1].set_xlabel("Classification", fontsize=12, fontweight='bold')
    fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical', fontsize=12, fontweight='bold')

    # Custom legend
    legend_handles = []
    for group_name, types in legend_labels.items():
        legend_handles.append(Patch(facecolor='white', edgecolor='white', label=f"**{group_name}**"))
        for t in types:
            legend_handles.append(Patch(facecolor=colors[t], label=f"  {t}"))

    fig.legend(
        handles=legend_handles,
        loc='center left',
        bbox_to_anchor=(0.85, 0.5),
        frameon=False
    )

    for text in fig.legends[0].get_texts():
        if text.get_text().startswith("**") and text.get_text().endswith("**"):
            text.set_text(text.get_text().strip("*"))
            text.set_weight('bold')
        else:
            text.set_fontstyle('normal')

    plt.tight_layout(rect=[0.07, 0, 0.85, 1])  # extra space on left for y-label and right for legend

    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def create_dist_plot(df, output_pdf):
    """
    Creates ordered vertical subplots of SV length distribution using plt.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        None: The plot is generated using plt.
    """
    # Log transform SV length
    df['log2_SV_len'] = np.log2(df['SVLEN'])

    # Define categories based on 'CLASSIFICATION'
    def categorize_sv_type(sv_type):
        if sv_type in ['HOMO']:
            return 'Tandem repeat'
        elif sv_type in ['STR']:
            return 'Tandem repeat'
        elif sv_type in ['TR']:
            return 'Tandem repeat'
        elif sv_type in ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']:
            return 'Mobile element'
        else:
            return 'NON_REPETITIVE'

    df['CLASSIFICATION_CATEGORY'] = df['CLASSIFICATION'].apply(categorize_sv_type)

    # Define the desired order of SV types
    ordered_classification_types = ['NON_REPETITIVE', 'HOMO', 'STR', 'TR', 'LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']

    # Define colors based on the provided image
    colors = {
        'NON_REPETITIVE': '#94d0c5',
        'HOMO': '#beb9d8',
        'STR': '#eb8274',
        'TR': '#87b0d2',
        'LTR': '#f4b567',
        'LINE': '#b4d56f',
        'SINE': '#f6cee0',
        'Retroposon': '#d9d9d8',
        'DNA': '#f9f6b7'
    }

    num_subplots = len(ordered_classification_types)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(15, 1 * num_subplots), sharex=True)
    fig.suptitle("SV length distribution for each classification type", fontsize=14, fontweight='bold')
    if num_subplots == 1:
        axes = [axes]

    legend_handles = []

    for i, classification_type in enumerate(ordered_classification_types):
        ax = axes[i]
        classification_type_df = df[df['CLASSIFICATION'] == classification_type]
        color = colors.get(classification_type, 'gray')

        if not classification_type_df.empty:
            x_min = 4
            x_max = df['log2_SV_len'].max()
            x_grid = np.linspace(x_min, x_max, 1000)

            if classification_type_df['CLASSIFICATION_CATEGORY'].iloc[0] == 'Mobile element':
                complete_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Full']
                fragment_df = classification_type_df[classification_type_df['RECIPROCAL'] == 'Partial']

                # KDE and normalization for "Full"
                if len(complete_df) > 1:
                    kde = stats.gaussian_kde(complete_df['log2_SV_len'])
                    y = kde(x_grid)
                    y /= y.max()  # Normalize peak to 1
                    ax.fill_between(x_grid, y, color=color, alpha=1, label=f'{classification_type} (Full)')
                    # Outline
                    ax.plot(x_grid, y, color='grey', linewidth=0.1)
                else:
                    print(f"Skipping {classification_type} (Full) due to insufficient data for KDE.")
                
                # KDE and normalization for "Partial"
                if len(fragment_df) > 1:
                    kde = stats.gaussian_kde(fragment_df['log2_SV_len'])
                    y = kde(x_grid)
                    y /= y.max()
                    ax.fill_between(x_grid, y, color=color, alpha=0.5, label=f'{classification_type} (Partial)')
                    ax.plot(x_grid, y, color='grey', linewidth=0.1)
                else:
                    print(f"Skipping {classification_type} (Partial) due to insufficient data for KDE.")

            else:
                kde = stats.gaussian_kde(classification_type_df['log2_SV_len'])
                y = kde(x_grid)
                y /= y.max()
                ax.fill_between(x_grid, y, color=color, alpha=1, label=classification_type)

            # Plot settings
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
            # xticks = [4, 8, 12, 16]
            # ax.set_xticks(xticks)
            # ax.xaxis.set_major_locator(plt.FixedLocator(xticks))
            # ax.set_xticklabels([str(x) for x in xticks])
            ax.set_xlabel('log$_2$(SV Length)', fontsize=12, fontweight='bold')
            ax.set_xticklabels([str(x) for x in custom_decimal_ticks])

            ax.xaxis.set_major_locator(plt.FixedLocator(custom_log2_ticks))

        for xtick in ax.get_xticks():
            ax.axvline(x=xtick, linestyle='--', color='gray', linewidth=0.5, zorder=1)
        
        # Add outline around plot
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_edgecolor('black')

    fig.text(0.1, 0.5, 'Density', va='center', rotation='vertical', fontsize=12, fontweight='bold')
    
    # Define SV type groups and colors for the legend
    legend_labels = {
        'NON_REPETITIVE': ['NON_REPETITIVE'],
        'Tandem repeat': ['HOMO', 'STR', 'TR'],
        'Mobile element': ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']
    }
    # Custom legend with bold section headers
    legend_handles = []
    for group_name, types in legend_labels.items():
        legend_handles.append(Patch(facecolor='white', edgecolor='white', label=f"**{group_name}**"))  # Group title
        for t in types:
            legend_handles.append(Patch(facecolor=colors[t], label=f"  {t}"))

    plt.subplots_adjust(right=0.75, hspace=0.9)  # Maintain the horizontal space

    legend = plt.legend(
        handles=legend_handles,
        loc='center left',
        bbox_to_anchor=(1.01, 8),  # Adjust the vertical position (y-coordinate)
        title=None,
        frameon=False
    )

    # Format bold section headers in the legend
    for text in legend.get_texts():
        if text.get_text().startswith("**") and text.get_text().endswith("**"):
            text.set_text(text.get_text().strip("*"))
            text.set_weight('bold')
        else:
            text.set_fontstyle('normal')

    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def merge_pdfs(output_pdf, input_pdfs):
    writer = PdfWriter()

    for pdf_file in input_pdfs:
        reader = PdfReader(pdf_file)
        for page in reader.pages:
            writer.add_page(page)

    # Save the merged PDF
    with open(output_pdf, "wb") as f_out:
        writer.write(f_out)

def main():
    # Read command-line arguments
    args = read_cmd_args()
    os.makedirs(args.out, exist_ok=True)

    # Read the TSV file
    df = read_tsv(args.tsv)

    output_pdf_0 = os.path.join(args.out, 'distributions.pdf')
    create_dist_plot(df.copy(), output_pdf_0)
    
    sv_order = ['INS', 'DEL', 'INV', 'DUP', 'BND']
    output_pdf_1 = os.path.join(args.out, 'histograms.pdf')
    create_hist_plot(df.copy(), output_pdf_1, sv_order)

    output_pdf = os.path.join(args.out, 'plots.pdf')
    merge_pdfs(output_pdf, [output_pdf_0, output_pdf_1])

if __name__ == "__main__":
    main()