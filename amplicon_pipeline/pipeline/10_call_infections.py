"""Call infections using chosen amplicons, coverage, and percent passing."""

import argparse
import os
from pathlib import Path

import pandas as pd

DEFAULT_AMPLICONS = ['ACS8_6_F_R', 'SERA2_F_R', 'AMA1_F_R']
DEFAULT_MIN_COVERAGE = 50
DEFAULT_PERCENT_AMPLICONS = 100

# input and output (exported by config.sh)
OUTPUT_DIR = Path(os.environ['OUTPUT_DIR'])
COVERAGE_TABLE = Path(os.environ['COVERAGE_TABLE'])
RESULTS_DIR = OUTPUT_DIR / 'infection_calling'
SAMPLE_LIST = RESULTS_DIR / 'filtered_samples.csv'


def load_analysis_coverage():
    coverage = pd.read_csv(COVERAGE_TABLE, sep='\t')
    samples = pd.read_csv(SAMPLE_LIST)

    coverage = coverage.rename(columns={'sample': 'seq_sample_id'})
    coverage['seq_sample_id'] = coverage['seq_sample_id'].astype(str)
    samples['seq_sample_id'] = samples['seq_sample_id'].astype(str)
    ids = set(samples['seq_sample_id'])
    coverage = coverage[coverage['seq_sample_id'].isin(ids)]
    return coverage


def call_infections(coverage_meta, amplicons, min_coverage, percent_amplicons):
    selected = coverage_meta[coverage_meta['amplicon'].isin(amplicons)]
    coverage = selected.pivot_table(
        index='seq_sample_id', columns='amplicon', values='mean_coverage'
    )
    coverage = coverage.reindex(
        index=sorted(coverage_meta['seq_sample_id'].unique()),
        columns=amplicons,
    )

    percent_passing = (coverage >= min_coverage).mean(axis=1) * 100
    infected = percent_passing >= percent_amplicons
    calls = coverage.reset_index()
    calls['infection_status'] = infected.map({True: 'INFECTED', False: 'NOT_INFECTED'}).values
    return calls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Call infection from the percentage of chosen amplicons above a coverage cutoff.'
    )
    parser.add_argument('--amplicons', nargs='+', default=DEFAULT_AMPLICONS,
                        metavar='AMPLICON', help='Amplicons used for calling.')
    parser.add_argument('--min-coverage', type=float, default=DEFAULT_MIN_COVERAGE,
                        help='Coverage required at each amplicon (default: 50).')
    parser.add_argument('--percent-amplicons', type=float, default=DEFAULT_PERCENT_AMPLICONS,
                        help='Percent of amplicons that must pass (default: 100).')
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    coverage_meta = load_analysis_coverage()
    infection_calls = call_infections(
        coverage_meta, args.amplicons, args.min_coverage, args.percent_amplicons
    )

    parameters = pd.DataFrame({
        'amplicon': args.amplicons,
        'min_coverage_required': args.min_coverage,
        'percent_amplicons_required': args.percent_amplicons,
    })
    parameters.to_csv(RESULTS_DIR / 'infection_calling_parameters.csv', index=False)
    infection_calls.to_csv(RESULTS_DIR / 'infection_calls.csv', index=False)

    print(f'Amplicons: {args.amplicons}')
    print(f'Coverage threshold: {args.min_coverage} reads')
    print(f'Amplicons required: {args.percent_amplicons}%')
    print(infection_calls['infection_status'].value_counts().to_string())
    print(f'Results written to: {RESULTS_DIR}')

