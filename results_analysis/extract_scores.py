import os
import re
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns


def parse_tm_score(file_path):
    our_tm, rmsd, ben_tm, num_subunits = None, None, None, None
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("Our best model:"):
                match = re.search(r'cb_(\d+)_output', line)
                if match:
                    num_subunits = int(match.group(1))
            elif line.startswith("TM-score:") and our_tm is None:
                our_tm = float(line.split(":")[1].strip())
            elif line.startswith("RMSD:"):
                rmsd = float(line.split(":")[1].strip())
            elif line.startswith("Ben's best model:"):
                next_line_idx = lines.index(line) + 1
                if next_line_idx < len(lines):
                    next_line = lines[next_line_idx]
                    if next_line.startswith("TM-score:"):
                        ben_tm = float(next_line.split(":")[1].strip())
    except Exception as e:
        print(f"Warning: Failed to parse {file_path} ({e})")
    return our_tm, rmsd, ben_tm, num_subunits

def main(main_dir):
    result = []
    subdirs = ['combfold', 'combfold_all', 'combfold_us_trivial', 'combfold_trivial']

    for complex_name in os.listdir(main_dir):
        complex_path = os.path.join(main_dir, complex_name)
        if not os.path.isdir(complex_path):
            continue

        row = {"Complex": complex_name}
        ben_tm_collected = None

        for subdir in subdirs:
            sub_path = os.path.join(complex_path, subdir)
            tm_file = os.path.join(sub_path, 'tm_score.txt')

            if not os.path.exists(tm_file):
                print(f"Warning: Missing file {tm_file}")
                continue

            our_tm, rmsd, ben_tm, num_subunits = parse_tm_score(tm_file)

            # Fill values
            row[f"{subdir}_TM"] = our_tm
            row[f"{subdir}_RMSD"] = rmsd
            row[f"{subdir}_num_subunits"] = num_subunits

            if ben_tm is not None:
                if ben_tm_collected is None:
                    ben_tm_collected = ben_tm
                elif abs(ben_tm - ben_tm_collected) > 1e-4:
                    print(f"Warning: Inconsistent Ben's TM-score in {tm_file}")

        row["Ben_best_TM"] = ben_tm_collected
        result.append(row)

    df = pd.DataFrame(result)
    df.to_csv(os.path.join(main_dir, "combined_scores.csv"), index=False)
    # ---- Plot: Bar Plot for combfold vs combfold_trivial ----
    label_map = {
        'combfold_TM': 'RedefineSubunit',
        'combfold_trivial_TM': 'CombFold'
    }
    bar_df = df[['combfold_TM', 'combfold_trivial_TM']].mean().reset_index()
    bar_df.columns = ['Method', 'Mean_TM']
    bar_df['Method'] = bar_df['Method'].map(label_map)

    plt.figure(figsize=(6, 4))
    sns.barplot(data=bar_df, x='Method', y='Mean_TM', palette=['#4C72B0', '#55A868'])
    plt.title('Mean TM-score: RedefineSubunit vs CombFold')
    plt.ylabel('TM-score')
    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, "mean_tm_score_barplot.png"))
    plt.close()

    # ---- Plot: Scatter Plot combfold vs combfold_trivial ----
    plt.figure(figsize=(6, 6))
    plt.scatter(df['combfold_trivial_TM'], df['combfold_TM'], color='#4C72B0')
    plt.plot([0, 1], [0, 1], '--', color='gray')
    plt.xlabel('CombFold TM-score')
    plt.ylabel('RedefineSubunit TM-score')
    plt.title('RedefineSubunit vs CombFold TM-score')
    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, "RedefineSubunit_vs_CombFold_scatter.png"))
    plt.close()

    # ---- Final Plot: Success Rate (High > 0.8, Acceptable > 0.7) ----
    success_data = []
    for method, label in [('combfold_trivial_TM', 'CF'), ('combfold_TM', 'ReDefineSubunit')]:
        scores = df[method].dropna()
        high = (scores > 0.8).sum()
        acceptable = (scores > 0.7).sum()
        total = len(scores)
        success_data.append({
            'Method': label,
            'High': high / total,
            'Acceptable': (acceptable - high) / total
        })

    success_df = pd.DataFrame(success_data)
    success_df.set_index('Method', inplace=True)

    # Define colors
    colors = {
        'CF': ['#D55E00', '#FFA07A'],  # Darker orange and lighter orange
        'ReDefineSubunit': ['#56B4E9', '#AEDFF7']  # Blue shades
    }

    # Plot
    # Plot
    fig, ax = plt.subplots(figsize=(2.5, 4))

    for idx, method in enumerate(success_df.index):
        high = success_df.loc[method, 'High']
        acceptable = success_df.loc[method, 'Acceptable']
        ax.bar(idx, high, color=colors[method][0], width=0.5)
        ax.bar(idx, acceptable, bottom=high, color=colors[method][1], width=0.5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['CF', 'ReDefineSubunit'], rotation=0)
    ax.set_ylabel('Success rate')
    ax.set_ylim(0, 1)
    ax.set_title('TM-score Success Rate', pad=10)
    ax.spines[['top', 'right']].set_visible(False)

    fig, ax = plt.subplots(figsize=(2.5, 4))

    for idx, method in enumerate(success_df.index):
        high = success_df.loc[method, 'High']
        acceptable = success_df.loc[method, 'Acceptable']
        ax.bar(idx, high, color=colors[method][0], width=0.5)
        ax.bar(idx, acceptable, bottom=high, color=colors[method][1], width=0.5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['CF', 'ReDefineSubunit'], rotation=0)
    ax.set_ylabel('Success rate')
    ax.set_ylim(0, 1)
    ax.set_title('TM-score Success Rate', pad=10)
    ax.spines[['top', 'right']].set_visible(False)

    # Add this block here (after plotting bars, before layout)
    fig, ax = plt.subplots(figsize=(2.5, 4))

    for idx, method in enumerate(success_df.index):
        high = success_df.loc[method, 'High']
        acceptable = success_df.loc[method, 'Acceptable']
        ax.bar(idx, high, color=colors[method][0], width=0.5)
        ax.bar(idx, acceptable, bottom=high, color=colors[method][1], width=0.5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['CF', 'ReDefineSubunit'], rotation=0)
    ax.set_ylabel('Success rate')
    ax.set_ylim(0, 1)
    ax.set_title('TM-score Success Rate', pad=10)
    ax.spines[['top', 'right']].set_visible(False)

    #  Add this block here (after plotting bars, before layout)
    # legend_elements = [
    #     plt.Rectangle((0, 0), 1, 1, color='#D55E00', label='High (CF)'),
    #     plt.Rectangle((0, 0), 1, 1, color='#FFA07A', label='Acceptable (CF)'),
    #     plt.Rectangle((0, 0), 1, 1, color='#0072B2', label='High (ReDefineSubunit)'),
    #     plt.Rectangle((0, 0), 1, 1, color='#56B4E9', label='Acceptable (ReDefineSubunit)'),
    # ]
    # ax.legend(handles=legend_elements)

    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, "tm_score_success_rate.png"))
    plt.close()

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, color='#D55E00', label='High (CF)'),
        plt.Rectangle((0, 0), 1, 1, color='#FFA07A', label='Acceptable (CF)'),
        plt.Rectangle((0, 0), 1, 1, color='#0072B2', label='High (RedefineSubunit)'),
        plt.Rectangle((0, 0), 1, 1, color='#56B4E9', label='Acceptable (RedefineSubunit)'),
    ]

    fig_legend = plt.figure(figsize=(3, 1.5))
    fig_legend.legend(handles=legend_elements, loc='center', frameon=False, ncol=2, borderaxespad=0.5)
    fig_legend.tight_layout(pad=1.0)
    fig_legend.savefig(os.path.join(main_dir, "tm_score_success_rate_legend.png"), dpi=300, bbox_inches='tight',
                       pad_inches=0.1)
    plt.close(fig_legend)

    print("Saved: combined_scores.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract TM-scores and RMSDs from CombFold results.")
    parser.add_argument("main_dir", help="Path to the main directory containing complex subdirectories.")
    args = parser.parse_args()
    main(args.main_dir)
