"""Filter samples for downstream infection analysis.

Inputs:
* Coverage table from 06_coveragetable.sh.
* Metadata
* Plate QC summary from 07_plateQC.py.
* Optional lists of samples to exclude.

Output in `infection_calling/`:
* `filtered_samples.csv` - sample IDs retained after filters (control, plate QC,
  HbS coverage, metadata, and duplicate-pairs).
"""

import os 
from pathlib import Path

import pandas as pd

# parameters 
HBS_AMPLICON = 'HBS_F_R'
MIN_HBS_COVERAGE = 10

# output
OUTPUT_DIR = Path(os.environ['OUTPUT_DIR'])

# inputs 
COVERAGE_TABLE = Path(os.environ['COVERAGE_TABLE'])
METADATA = Path(os.environ['METADATA'])
# optional filter lists: will skip if file is not there
NO_METADATA_LIST = os.environ.get('NO_METADATA_LIST', '')
DUPLICATES_LIST = os.environ.get('DUPLICATES_LIST', '')
PLATE_QC = OUTPUT_DIR / 'plateQC' / 'plate_qc_summary.csv'
RESULTS_DIR = OUTPUT_DIR / 'infection_calling'


def report_filter(step, description, before, after):
    print(f'Step {step} - {description}: {after}/{before} samples kept ({before - after} removed)')


def load_optional_list(path):
    """Read an optional CSV filter list; return None if unavailable or empty."""
    if not Path(path).is_file():
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None


if __name__ == '__main__':
    coverage = pd.read_csv(COVERAGE_TABLE, sep='\t')
    metadata = pd.read_csv(METADATA)
    plate_qc = pd.read_csv(PLATE_QC)
    no_metadata = load_optional_list(NO_METADATA_LIST)
    duplicates = load_optional_list(DUPLICATES_LIST)

    coverage = coverage.rename(columns={'sample': 'seq_sample_id'})
    coverage['seq_sample_id'] = coverage['seq_sample_id'].astype(str)
    metadata['seq_sample_id'] = metadata['seq_sample_id'].astype(str)
    metadata['plateID'] = metadata['plateID'].astype(str)
    plate_qc['plateID'] = plate_qc['plateID'].astype(str)

    ids = set(metadata['seq_sample_id'])

    before = len(ids)
    ids &= set(metadata.loc[metadata['sample_type'] == 'sample', 'seq_sample_id'])
    report_filter(1, 'Remove controls', before, len(ids))


    before = len(ids)
    if no_metadata is not None:
        ids -= set(no_metadata['seq_sample_id'].astype(str))
    report_filter(2, 'Remove samples lacking metadata', before, len(ids))


    before = len(ids)
    excluded_plates = set(plate_qc.loc[plate_qc['plate_status'] == 'EXCLUDE', 'plateID'])
    ids &= set(metadata.loc[~metadata['plateID'].isin(excluded_plates), 'seq_sample_id'])
    report_filter(3, 'Remove QC-failed plates', before, len(ids))


    hbs = coverage[coverage['amplicon'] == HBS_AMPLICON]
    if hbs.empty:
        raise ValueError(f"No coverage rows found for HBS amplicon '{HBS_AMPLICON}'")
    before = len(ids)
    ids &= set(hbs.loc[hbs['mean_coverage'] >= MIN_HBS_COVERAGE, 'seq_sample_id'])
    report_filter(4, f'Remove HBS coverage < {MIN_HBS_COVERAGE}', before, len(ids))


    duplicate_ids = set()
    if duplicates is not None:
        for _, row in duplicates.iterrows():
            id_1 = str(row['seq_sample_id_1'])
            id_2 = str(row['seq_sample_id_2'])
            if id_1 in ids and id_2 in ids:
                duplicate_ids.update([id_1, id_2])
    before = len(ids)
    ids -= duplicate_ids
    report_filter(5, 'Remove duplicate sample pairs', before, len(ids))


    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RESULTS_DIR / 'filtered_samples.csv'
    pd.DataFrame({'seq_sample_id': sorted(ids)}).to_csv(output_path, index=False)
    print(f'QC-passed sample list written to: {output_path}')
