# Amplicon sequencing pipeline

Pipeline for processing multiplexed amplicon sequencing data targeting human and
*Plasmodium falciparum* (Pf) loci, applied to 2,246 samples from Cameroon in
[Hopson et al. (2025)](https://doi.org/10.1101/2025.09.09.675015).

**Workflow**

1. **Process, align, and count** (scripts 01-06): FASTQ files to coverage depth summaries.
2. **QC and infection calling** (scripts 07-10): use coverage depth to flag plates, filter samples, and call infections.
3. **Variant calling** (scripts 11-14): perform joint variant calling across human and Pf amplicons.

## Setup

### Requirements

- Slurm, conda (conda must be on your PATH)
- Tools (cutadapt, FastQC, MultiQC, bwa, samtools, bcftools) and Python packages (pandas, matplotlib, scipy), all specified in `environment.yml`

### Reference genome

The pipeline aligns to a single FASTA containing both genomes:

- Human: [GRCh38 primary assembly, Ensembl release 108](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
- Parasite: [Plasmodium falciparum 3D7, PlasmoDB release 62](https://plasmodb.org/a/service/raw-files/release-62/Pfalciparum3D7/fasta/data/PlasmoDB-62_Pfalciparum3D7_Genome.fasta)

Concatenate the two FASTA files before running the pipeline, and set `REF` in
`config.sh` to the result. Contig names must match those in `data/bed/`.

### Inputs

Illumina MiSeq paired-end FASTQ files, one folder per sequencing run (one 384-sample plate in this study), named `<run>R`, with FASTQ files named `<run>X<number>`:

```
data/runs/24118R/Fastq/24118X100_..._R{1,2}_001.fastq.gz
```

Reads for this study are deposited in the NCBI SRA under BioProject
[PRJNA1503656](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1503656).

Copy `data/sample_metadata_template.csv` to `data/sample_metadata.csv` and add one row per sample.
Columns:

- `seq_sample_id` - FASTQ file name up to the first `_` (e.g. 24118X100)
- `plateID` - sequencing run (e.g. 24118)
- `sample_type` - `sample` or `control`
- `control_type` - `negative` or `positive`
- `Parasitemia` -  microscopy category (e.g. 0, Pf+); optional, only used in`09_select_infection_parameters.py`

### Configuration

Edit `config.sh`:

- `REF` - path to the concatenated reference (see above)
- `ALL_RUNS` - names of the run folders inside `data/runs/`, e.g. 24118R
- `SLURM_ACCOUNT`, `SLURM_PARTITION` - cluster account and partition
- `DATA_DIR` - FASTQ files, in subfolders by sequencing run (default `data/runs/`)

## Quick start

Run from `amplicon_pipeline/` after editing `config.sh`.

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate amplicon-pipeline
source ./config.sh
```

Create the BWA and samtools indexes:

```bash
bwa index "$REF"
samtools faidx "$REF"
```

### 2. Process FASTQ files and calculate coverage

```bash
bash submit_alignment.sh
```

Wait for this step to finish and create `output/all_runs_mean_coverage.txt`.

### 3. Perform plate QC

Summarize and plot coverage depth across controls:

```bash
python pipeline/07_plateQC.py
```

Select a cutoff based on plots, then filter samples:

```bash
python pipeline/07_plateQC.py --cutoff 0.02
python pipeline/08_filter_samples.py
```

### 4. Select parameters for infection calling

```bash
python pipeline/09_select_infection_parameters.py
```

Inspect the outputs and select the amplicons, minimum coverage depth, and
percentage of passing amplicons to use in the next step.

### 5. Call infections

```bash
python pipeline/10_call_infections.py \
    --amplicons ACS8_6_F_R SERA2_F_R AMA1_F_R --min-coverage 50 \
    --percent-amplicons 100
```

### 6. Call human and parasite variants

```bash
bash submit_variant_calling.sh
```

## Repository layout

```
config.sh                     pipeline settings
environment.yml
submit_alignment.sh           submits scripts 01–06 to Slurm with job dependencies
submit_variant_calling.sh     submits scripts 11–14 to Slurm; run after scripts 07–10
pipeline/                     numbered scripts, one per step
data/
  bed/                        human_pf_amplicons.bed (all amplicons) and one <amplicon>.bed per amplicon
  sample_metadata_template.csv   template for sample metadata
  sample_metadata.csv         sample metadata (see Inputs)
  runs/                       FASTQ files, one folder per sequencing run (see Inputs)
output/
  runs/                        trimmed reads, FastQC, BAMs, pileups, one folder per run
  all_runs_mean_coverage.txt   coverage depth per amplicon/sample
  multiqc/                     MultiQC report
  plateQC/                     control ratios, plots, plate pass/fail
  infection_calling/           filtered samples, parameter tables, infection calls
  all_bams.txt, infected_bams.txt   BAM lists for variant calling
  human/                       filtered VCFs and genotype tables
  parasite/                    filtered VCFs and genotype tables
logs/                          stdout/stderr
```


## Pipeline steps

### FASTQ processing and coverage: scripts 01–06

Run with `bash submit_alignment.sh`.

- `01_trim.sh` - trims adapters with cutadapt. Output: `output/runs/<run>/trimmed.fastq/`
- `02_fastqc.sh` - assesses trimmed-read quality. Output: `output/runs/<run>/fastqc/`
- `02b_multiqc.sh` - combines FastQC reports. Output: `output/multiqc/`
- `03_align.sh` - aligns reads to the joint human–Pf reference. Output: `output/runs/<run>/bams/`
- `04_pileup.sh` - calculates depth across amplicons. Output: `output/runs/<run>/mpileup_depth/<sample>.pileup.txt`
- `05_meancoverage.sh` - calculates mean coverage per sample and amplicon. Output: `output/runs/<run>/mpileup_depth/<sample>.mean_coverage.txt`
- `06_coveragetable.sh` - combines coverage across sequencing runs. Output: `output/all_runs_mean_coverage.txt`

### QC and infection calling: scripts 07–10

- `07_plateQC.py` - calculates negative-to-positive control coverage ratios.
  Run it first without a cutoff to inspect the results, then again with the
  selected cutoff to classify plates.
- `08_filter_samples.py` - removes controls, QC-failed plates, samples with HbS
  coverage depth below 10, samples lacking metadata, and duplicate pairs. The missing metadata and duplicate exclusion lists are optional.
- `09_select_infection_parameters.py` - compares parasite amplicons and coverage
  thresholds using microscopy parasitemia.
- `10_call_infections.py` - calls a sample infected when the specified
  percentage of selected amplicons meets the minimum coverage depth threshold.

Outputs are written to `output/plateQC/` and
`output/infection_calling/`.

### Variant calling: scripts 11–14

Run with `bash submit_variant_calling.sh`.

- `11_list_bams.sh` - creates BAM lists for all QC-passed samples and infected
  samples.
- `12_call_human_variants.sh` - jointly calls and filters human variants across
  QC-passed samples.
- `13_call_parasite_variants.sh` - jointly calls and filters parasite variants
  across infected samples.
- `14_format_vcf.sh` - optionally removes febrile samples, recomputes INFO tags,
  and combines the parasite amplicon VCFs.

Final VCFs are written to `output/human/<amplicon>.human.filtered.final.vcf.gz` and
`output/parasite/pf.all.filtered.final.vcf.gz`. Each VCF has a corresponding `.gts.txt` genotype table.
