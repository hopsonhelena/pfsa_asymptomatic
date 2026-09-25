"""Compare candidate amplicons and coverage depth cutoffs."""

import argparse
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

CANDIDATE_AMPLICONS = ['1127000_F_R', 'ACS82_2_F_2_R', 'ACS8_6_F_R',
                        'AMA1_F_R', 'CSP_F_R', 'SERA2_F_R', 'TRAP_F_R']
N_AMPLICONS_TO_RECOMMEND = 3
DEFAULT_REFERENCE_DEPTH = 50
THRESHOLD_PARASITEMIA_LEVEL = 'Pf+'
CORRELATION_SEPARATE_RANKS = 4
COVERAGE_DEPTH_THRESHOLDS = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]

# input and output (exported by config.sh)
OUTPUT_DIR = Path(os.environ['OUTPUT_DIR'])
COVERAGE_TABLE = Path(os.environ['COVERAGE_TABLE'])
METADATA = Path(os.environ['METADATA'])
RESULTS_DIR = OUTPUT_DIR / 'infection_calling'
SAMPLE_LIST = RESULTS_DIR / 'filtered_samples.csv'

def normalize_parasitemia(value):
    if isinstance(value, str):
        value = re.sub(r'[\s\u200b\u200c\u200d\ufeff]+', '', value)
        if value == '0':
            return 0
    return value


def load_analysis_data():
    coverage = pd.read_csv(COVERAGE_TABLE, sep='\t')
    metadata = pd.read_csv(METADATA)
    samples = pd.read_csv(SAMPLE_LIST)

    coverage = coverage.rename(columns={'sample': 'seq_sample_id'})
    coverage['seq_sample_id'] = coverage['seq_sample_id'].astype(str)
    metadata['seq_sample_id'] = metadata['seq_sample_id'].astype(str)
    metadata['Parasitemia'] = metadata['Parasitemia'].apply(normalize_parasitemia)
    samples['seq_sample_id'] = samples['seq_sample_id'].astype(str)
    ids = set(samples['seq_sample_id'])

    combined = coverage.merge(metadata, on='seq_sample_id', how='inner')
    combined = combined[combined['seq_sample_id'].isin(ids)]
    return combined


def infer_parasitemia_order(values):
    present = pd.Series(values).dropna().unique()

    def sort_key(value):
        if value == 0 or value == '0':
            return 0
        if isinstance(value, str) and value.startswith('Pf'):
            return value.count('+')
        return None

    keyed = [(sort_key(value), value) for value in present]
    keyed = [(key, value) for key, value in keyed if key is not None]
    keyed.sort(key=lambda pair: pair[0])

    return [value for key, value in keyed]


def compute_pf_plus_medians(coverage_meta):
    candidates = coverage_meta[coverage_meta['amplicon'].isin(CANDIDATE_AMPLICONS)]
    at_pf_plus = candidates[candidates['Parasitemia'] == THRESHOLD_PARASITEMIA_LEVEL]
    return at_pf_plus.groupby('amplicon')['mean_coverage'].median()


def select_amplicons(coverage_meta, threshold):
    candidates = coverage_meta[coverage_meta['amplicon'].isin(CANDIDATE_AMPLICONS)]
    calibration = candidates[candidates['Parasitemia'] == THRESHOLD_PARASITEMIA_LEVEL]
    prop_above_threshold = (
        calibration.groupby('amplicon')['mean_coverage']
        .apply(lambda x: (x >= threshold).mean())
        .sort_values(ascending=False)
    )

    recommended = prop_above_threshold.head(N_AMPLICONS_TO_RECOMMEND).index.tolist()
    return recommended, prop_above_threshold


def compute_coverage_by_group(coverage_meta, parasitemia_order):
    candidates = coverage_meta[coverage_meta['amplicon'].isin(CANDIDATE_AMPLICONS)]
    candidates = candidates[candidates['Parasitemia'].isin(parasitemia_order)].copy()

    level_rank = build_grouped_parasitemia_rank(parasitemia_order)
    candidates['parasitemia_level'] = candidates['Parasitemia'].map(level_rank)

    summary = (
        candidates.groupby(['amplicon', 'parasitemia_level'])['mean_coverage']
        .agg(['mean', 'median', 'count'])
        .reset_index()
    )
    return summary.sort_values(['amplicon', 'parasitemia_level']).reset_index(drop=True)


def build_grouped_parasitemia_rank(parasitemia_order, separate_ranks=CORRELATION_SEPARATE_RANKS):
    rank = {}
    for i, level in enumerate(parasitemia_order):
        rank[level] = min(i, separate_ranks)
    return rank


def compute_amplicon_correlation(coverage_meta, parasitemia_order):
    parasitemia_rank = build_grouped_parasitemia_rank(parasitemia_order)

    candidates = coverage_meta[coverage_meta['amplicon'].isin(CANDIDATE_AMPLICONS)]
    candidates = candidates[candidates['Parasitemia'].isin(parasitemia_order)].copy()
    candidates['parasitemia_rank'] = candidates['Parasitemia'].map(parasitemia_rank)

    rows = []
    for amplicon, group in candidates.groupby('amplicon'):
        rho, p_value = spearmanr(group['parasitemia_rank'], group['mean_coverage'])
        rows.append({'amplicon': amplicon, 'spearman_rho': rho, 'p_value': p_value, 'n': len(group)})

    return pd.DataFrame(rows).sort_values('spearman_rho', ascending=False).reset_index(drop=True)


def compute_amplicon_coverage_scale(coverage_meta, amplicons=CANDIDATE_AMPLICONS):
    candidates = coverage_meta[coverage_meta['amplicon'].isin(amplicons)]

    scale = candidates.groupby('amplicon')['mean_coverage'].agg(['mean', 'median', 'std', 'count'])
    scale['cv'] = scale['std'] / scale['mean']
    return scale.reset_index().sort_values('cv')


def compute_coverage_depth_curve(coverage_meta, amplicons, parasitemia_level=THRESHOLD_PARASITEMIA_LEVEL,
                                  thresholds=COVERAGE_DEPTH_THRESHOLDS):
    columns = ['amplicon', 'min_depth', 'proportion_above', 'n']
    subset = coverage_meta[
        coverage_meta['amplicon'].isin(amplicons) & (coverage_meta['Parasitemia'] == parasitemia_level)
    ]

    rows = []
    for amplicon, group in subset.groupby('amplicon'):
        for min_depth in thresholds:
            rows.append({
                'amplicon': amplicon,
                'min_depth': min_depth,
                'proportion_above': (group['mean_coverage'] >= min_depth).mean(),
                'n': len(group),
            })

    return pd.DataFrame(rows, columns=columns)


def plot_amplicon_ranking_elbow(amplicon_ranking, n_selected, output_path):
    ranked = amplicon_ranking.sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(range(1, len(ranked) + 1), ranked.values, marker='o')
    ax.axvline(n_selected + 0.5, color='red', linestyle='--', linewidth=1,
               label=f'cutoff after {n_selected} amplicons')

    ax.set_xticks(range(1, len(ranked) + 1))
    ax.set_xticklabels(ranked.index, rotation=45, ha='right')
    ax.set_xlabel('Amplicon (ranked)')
    ax.set_ylabel('Proportion of Pf+ samples above threshold')
    ax.set_title('Amplicon selection ranking')
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

def plot_coverage_depth_curve(depth_curve, highlighted_amplicons, output_path,
                              reference_depths=(DEFAULT_REFERENCE_DEPTH,)):
    fig, ax = plt.subplots(figsize=(8, 6))

    for amplicon in sorted(depth_curve['amplicon'].unique()):
        amp_data = depth_curve[depth_curve['amplicon'] == amplicon].sort_values('min_depth')

        if amplicon in highlighted_amplicons:
            ax.plot(amp_data['min_depth'], amp_data['proportion_above'], marker='o', markersize=4,
                     linewidth=2.5, label=f'{amplicon} (recommended)')
        else:
            ax.plot(amp_data['min_depth'], amp_data['proportion_above'], marker='o', markersize=3,
                     linewidth=1, label=amplicon)

    for reference_depth in reference_depths:
        ax.axvline(reference_depth, color='red', linestyle='--', linewidth=1)

    ax.set_xlabel('Minimum coverage depth (reads)')
    ax.set_ylabel('Proportion of Pf+ samples above depth')
    ax.set_title('Coverage depth curve (recommended amplicons highlighted)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Evaluate candidate amplicons and coverage thresholds.'
    )
    parser.add_argument('--reference-depth', type=float, default=DEFAULT_REFERENCE_DEPTH,
                        help='Coverage depth used to rank candidates and mark plots (default: 50).')
    args = parser.parse_args()
    reference_depth = args.reference_depth

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    coverage_meta = load_analysis_data()

    median_at_pf_plus = compute_pf_plus_medians(coverage_meta)
    recommended_amplicons, amplicon_ranking = select_amplicons(coverage_meta, reference_depth)

    parasitemia_order = infer_parasitemia_order(coverage_meta['Parasitemia'])
    print(f'Parasitemia order: {parasitemia_order}')
    print(f'Reference coverage depth: {reference_depth}')
    print(f'Top three amplicons: {recommended_amplicons}')
    print()
    print(f'Median coverage at {THRESHOLD_PARASITEMIA_LEVEL}:')
    print(median_at_pf_plus.sort_values(ascending=False).to_string())
    print()
    print('Proportion of Pf+ samples above threshold:')
    print(amplicon_ranking.to_string())

    median_at_pf_plus.sort_values(ascending=False).rename('median_coverage').to_csv(
        RESULTS_DIR / 'median_coverage_at_pf_plus.csv'
    )
    amplicon_ranking.rename('prop_above_threshold').to_csv(
        RESULTS_DIR / 'amplicon_selection_ranking.csv'
    )

    coverage_by_group = compute_coverage_by_group(coverage_meta, parasitemia_order)
    coverage_by_group.to_csv(RESULTS_DIR / 'coverage_by_parasitemia_group.csv', index=False)

    correlation = compute_amplicon_correlation(coverage_meta, parasitemia_order)
    correlation.to_csv(RESULTS_DIR / 'amplicon_parasitemia_correlation.csv', index=False)

    coverage_scale = compute_amplicon_coverage_scale(coverage_meta)
    coverage_scale.to_csv(RESULTS_DIR / 'amplicon_coverage_scale.csv', index=False)

    depth_curve = compute_coverage_depth_curve(coverage_meta, CANDIDATE_AMPLICONS)
    depth_curve.to_csv(RESULTS_DIR / 'coverage_depth_curve_all_candidates.csv', index=False)

    ranking_elbow_plot_path = RESULTS_DIR / 'amplicon_ranking_elbow.png'
    plot_amplicon_ranking_elbow(amplicon_ranking, N_AMPLICONS_TO_RECOMMEND, ranking_elbow_plot_path)

    depth_curve_plot_path = RESULTS_DIR / 'coverage_depth_curve.png'
    plot_coverage_depth_curve(depth_curve, recommended_amplicons, depth_curve_plot_path,
                              reference_depths=(reference_depth,))

    print()
    print('Coverage by amplicon (lower CV is more consistent):')
    print(coverage_scale.to_string(index=False))
    print()
    print('Spearman correlation with parasitemia rank:')
    print(correlation.to_string(index=False))
    print()
    print(f'Proportion of Pf+ samples above reference depth = {reference_depth}:')
    retention_at_reference = depth_curve[depth_curve['min_depth'] == reference_depth]
    retention_at_reference = retention_at_reference.sort_values('proportion_above', ascending=False)
    print(retention_at_reference[['amplicon', 'proportion_above']].to_string(index=False))

    print(f'\nResults written to: {RESULTS_DIR}')

