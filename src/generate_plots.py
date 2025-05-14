import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def read_cmd_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot SV types by repeat category and SV length distribution.")
    parser.add_argument('--tsv', required=True, help='Path to input TSV file')
    parser.add_argument('--out', required=True, help='Output dir')
    return parser.parse_args()

def read_tsv(tsv_file):
    """Read and preprocess the TSV file."""
    df = pd.read_csv(tsv_file, sep='\t')
    df['SV_len'] = pd.to_numeric(df['SV_len'], errors='coerce')
    df = df.dropna(subset=['SV_len'])
    return df

def create_hist_plot(df, output_dir):
    """Generate the plot for SV types.

    Returns:
        None: The plot is generated using plt.
    """
    # Define SV type groups and colors
    legend_labels = {
        'Non-repetitive': ['Non-repetitive'],
        'Tandem repeat': ['HOMO', 'STR', 'TR'],
        'Mobile element': ['LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']
    }

    # Flatten the hue order
    svtype_order = [svtype for group in legend_labels.values() for svtype in group]

    # Set color map
    colors = {
        'Non-repetitive': '#94d0c5',
        'HOMO': '#beb9d8',
        'STR': '#eb8274',
        'TR': '#87b0d2',
        'LTR': '#f4b567',
        'LINE': '#b4d56f',
        'SINE': '#f6cee0',
        'Retroposon': '#d9d9d8',
        'DNA': '#f9f6b7'
    }

    # Group and count
    grouped = df.groupby(['SV', 'SV_type']).size().reset_index(name='Count')

    # Enforce SV type and SV class order
    grouped['SV_type'] = pd.Categorical(grouped['SV_type'], categories=svtype_order, ordered=True)
    sv_order = ['INS', 'DEL', 'INV', 'DUP', 'BND']
    grouped['SV'] = pd.Categorical(grouped['SV'], categories=sv_order, ordered=True)

    # Plot
    fig = plt.figure(figsize=(10, 6))
    sns.set(style="whitegrid")
    ax = sns.barplot(data=grouped, x='SV', y='Count', hue='SV_type', hue_order=svtype_order, palette=colors)

    # Custom legend with bold section headers
    legend_handles = []
    for group_name, types in legend_labels.items():
        legend_handles.append(Patch(facecolor='white', edgecolor='white', label=f"**{group_name}**"))  # Group title
        for t in types:
            legend_handles.append(Patch(facecolor=colors[t], label=f"  {t}"))

    legend = ax.legend(
        handles=legend_handles,
        loc='center left',
        bbox_to_anchor=(1.0, 0.5),
        title=None,
        frameon=False
    )

    # Format bold section headers
    for text in legend.get_texts():
        if text.get_text().startswith("**") and text.get_text().endswith("**"):
            text.set_text(text.get_text().strip("*"))
            text.set_weight('bold')
        else:
            text.set_fontstyle('normal')

    # Axis labels and layout
    ax.set_title("Structural Variants by Repeat Category")
    ax.set_xlabel("SV Class")
    ax.set_ylabel("Count")
    plt.tight_layout()

    output_pdf = os.path.join(output_dir, 'histograms.pdf')
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def create_dist_plot(df, output_dir):
    """
    Creates ordered vertical subplots of SV length distribution using plt.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        None: The plot is generated using plt.
    """
    # Map STATUS to Complete and Fragment
    # df['STATUS_RECODED'] = df['STATUS'].replace({'Full': 'Complete', 'Partial': 'Fragment', 'Minimal': 'Fragment'})

    # Log transform SV length
    df['log2_SV_len'] = np.log2(df['SV_len'])

    # Define categories based on 'SV_type'
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
            return 'Non-repetitive'

    df['SV_CATEGORY'] = df['SV_type'].apply(categorize_sv_type)

    # Define the desired order of SV types
    ordered_sv_types = ['Non-repetitive', 'HOMO', 'STR', 'TR', 'LTR', 'LINE', 'SINE', 'Retroposon', 'DNA']

    # Define colors based on the provided image
    colors = {
        'Non-repetitive': '#94d0c5',
        'HOMO': '#beb9d8',
        'STR': '#eb8274',
        'TR': '#87b0d2',
        'LTR': '#f4b567',
        'LINE': '#b4d56f',
        'SINE': '#f6cee0',
        'Retroposon': '#d9d9d8',
        'DNA': '#f9f6b7'
    }

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    num_subplots = len(ordered_sv_types)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(5, 0.8 * num_subplots), sharex=True)
    if num_subplots == 1:
        axes = [axes]

    legend_handles = []

    for i, sv_type in enumerate(ordered_sv_types):
        ax = axes[i]
        sv_type_df = df[df['SV_type'] == sv_type]
        color = colors.get(sv_type, 'gray')

        if not sv_type_df.empty:
            if sv_type_df['SV_CATEGORY'].iloc[0] == 'Mobile element':
                complete_df = sv_type_df[sv_type_df['STATUS'] == 'Full']
                fragment_df = sv_type_df[sv_type_df['STATUS'] == 'Partial']
                sns.kdeplot(data=complete_df, x='log2_SV_len', fill=True, alpha=1, color=color, label=f'{sv_type} (Full)', ax=ax)
                sns.kdeplot(data=fragment_df, x='log2_SV_len', fill=True, alpha=0.3, color=color, label=f'{sv_type} (Partial)', ax=ax)
            else:
                sns.kdeplot(data=sv_type_df, x='log2_SV_len', fill=True, alpha=1, color=color, label=sv_type, ax=ax)

            ax.set_yticks([0, 1])
            ax.set_yticklabels(['0', '1'])
            ax.set_ylim(0, 1)
            ax.set_xlim(4, df['log2_SV_len'].max())
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.tick_params(axis='y', which='major', length=5)
            ax.tick_params(axis='y', which='minor', length=2)
            ax.yaxis.set_minor_locator(plt.FixedLocator([0.5]))
            ax.set_xticks([4, 8, 12, 16])
            ax.set_ylabel('') # Ensure no y-axis label is set

        else:
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['0', '1'])
            ax.set_ylim(0, 1)
            ax.set_xlim(4, df['log2_SV_len'].max())
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.set_xticks([4, 8, 12, 16])
            ax.tick_params(axis='y', which='major', length=5)
            ax.tick_params(axis='y', which='minor', length=2)
            ax.yaxis.set_minor_locator(plt.FixedLocator([0.5]))
            ax.spines['bottom'].set_visible(False)
            # ax.set_ylabel('') # Ensure no y-axis label is set

        if i < num_subplots - 1:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('log$_2$(SV Length)')

    # Define SV type groups and colors for the legend
    legend_labels = {
        'Non-repetitive': ['Non-repetitive'],
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


    output_pdf = os.path.join(output_dir, 'distributions.pdf')
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def main():
    # Read command-line arguments
    args = read_cmd_args()
    os.makedirs(args.out, exist_ok=True)

    # Read the TSV file
    df = read_tsv(args.tsv)

    # Generate the plots using plt
    create_dist_plot(df.copy(), args.out)
    create_hist_plot(df.copy(), args.out)

if __name__ == "__main__":
    main()