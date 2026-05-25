# OctoPan — Pangenome-Based WES/WGS Analysis Pipeline

**OctoPan** is a comprehensive, scalable WDL (Workflow Description Language) pipeline for whole-exome sequencing (WES) and whole-genome sequencing (WGS) analysis. It leverages a **personalised human pangenome reference graph** (HPRC v1.1) to maximise read mapping sensitivity, then surjects alignments back to linear references (GRCh38 and CHM13) for downstream variant calling, structural variant detection, STR genotyping, mitochondrial analysis, and splicing-variant prioritisation — all in a single, automated run.

---
<img width="1580" height="887" alt="image" src="https://github.com/user-attachments/assets/1eec3743-29e3-4304-93bd-88b65a13749a" />

---
## Why OctoPan?

Traditional pipelines align reads to a single linear reference (GRCh38), which introduces reference bias — reads carrying alternative alleles or from underrepresented populations are less likely to map correctly. OctoPan addresses this by:

- **Mapping to a pangenome graph** built from 94 diverse haplotypes (HPRC v1.1), substantially reducing reference bias and improving variant recall in complex and repetitive regions.
- **Dual-reference surjection** — aligned reads are projected onto both GRCh38 and CHM13 linear references, enabling variant calling with two independent assemblies for cross-validation.
- **Integrated end-to-end analysis** — a single pipeline covers QC, alignment, SNV/indel calling (DeepVariant), structural variant calling (Manta), short-tandem repeat expansion analysis (STRipy), mitochondrial variant analysis (MitoHPC), splicing variant prioritisation (SQUIRLS), mutational signature analysis (SNVstory), and sex determination.
- **Containerised tools** — all bioinformatics tools are packaged in Singularity images, so no manual installation of complex software is required.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Repository Structure](#2-repository-structure)
3. [Input Files](#3-input-files)
4. [Configuring the JSON Input File](#4-configuring-the-json-input-file)
5. [Running the Pipeline](#5-running-the-pipeline)
6. [Workflow Steps](#6-workflow-steps)
7. [Outputs](#7-outputs)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

### Java (Required)

OctoPan runs on [Cromwell](https://cromwell.readthedocs.io/), which requires Java. Ensure the following version (or equivalent OpenJDK 18+) is available in your environment:

```bash
java -version
```

### R Dependencies

The STRipy post-processing script (`json_to_csv_graph.R`) requires the following R packages:

```r
install.packages(c("jsonlite", "purrr", "dplyr"))
```

### Python Dependencies

The sex-determination script (`determine_sex_from_bam.py`) requires:

```bash
pip install pandas
```

### Reqiered Datasets

#### genomAD
you need to get `gnomad.exomes.v4.0.sites.chr1-22,X,Y.vcf.gz`. Firstly, download all exomes sites in VCF formats, from chr1 to chr22, X, and Y from genomAD portal here
```
https://gnomad.broadinstitute.org/data#v4-variants
```
Then, use `bcftools concat` to combine all of them into one VCF file and save it with this name `gnomad.exomes.v4.0.sites.chr1-22,X,Y.vcf.gz`. Finaly, place the VCF inside the `datasets` directory.

#### snvstory_resources
download it from Zenodo, unzip and place it inside the `containers` directory

```
https://zenodo.org/records/20378382/files/snvstory_resource.tar.gz?download=1
```

#### Other datasets
Other required datasets are stored on Zenodo. Just download them, unzip, and store them inside the `datasets` directory.

```
wget https://zenodo.org/records/20378382/files/datasets.tar.gz?download=1
```

### Singularity

All other tools (DeepVariant, STRipy, MitoHPC, SQUIRLS, SNVstory, OctoPan toolkit) are distributed as Singularity images stored on Zenodo.
download them from this link 

```
https://zenodo.org/records/20378382/files/singularities.tar.gz?download=1
```
Unpack them and store them inside the `containers/` directory. No further installation is needed for these tools.

---

## 2. Repository Structure
The directories structure after downloading all the dependencies should follow this:

```
OctoPan/
├── Main.wdl                  # Top-level workflow entry point
├── Main_inputs.json          # JSON file with all input parameters (edit this)
├── WDL_v3.1_fixed.conf       # Cromwell backend configuration
├── cromwell_slurm.conf       # SLURM-specific Cromwell configuration
├── cromwell-87.jar           # Cromwell execution engine
│
├── wdls/                     # All sub-workflow WDL modules
│   ├── OctoPan.wdl           # Main orchestrating sub-workflow
│   ├── QualityControl.wdl
│   ├── MapToPangenome.wdl
│   ├── SurjectBAM.wdl
│   ├── VariantCalling.wdl
│   ├── Manta.wdl
│   ├── STRipyPipeline.wdl
│   ├── SQUIRLS.wdl
│   ├── SNVstory.wdl
│   ├── ExtractChrM.wdl
│   ├── SexDetermination.wdl
│   ├── MitoHPC.wdl
│   ├── IndexKmers.wdl
│   └── Versions.wdl
│
├── containers/               # Singularity images and helper scripts
│   ├── octopan-toolkit.sif
│   ├── octopan-toolkit.def
│   ├── google_deepvariant_1.6.0.sif
│   ├── stripy/
│   ├── stripy.sif
│   ├── mitohpc.sif
│   ├── snvstory.sif
│   ├── snvstory_resource/
│   ├── squirls/
│   ├── snpEff/
│   ├── determine_sex_from_bam.py
│   └── json_to_csv_graph.R
│
└── datasets/                 # Reference files and annotation databases
    ├── hprc-v1.1-mc-chm13.gbz
    ├── hprc-v1.1-mc-chm13.hapl
    ├── grch38.paths.txt
    ├── chm13.paths.txt
    ├── GRCh38.dict
    ├── GRCh38_edit_chr_names.fa
    ├── CHM13_edit_chr_names.fa
    ├── Homo_sapiens_assembly38.fasta
    ├── disease_superpanel_example.bed
    └── gnomad.exomes.v4.0.sites.chr1-22,X,Y.vcf.gz
```
you've to make sure that you downloaded all datasets and containers and placed them correctly 

---

## 3. Input Files

### 3.1 FASTQ Files

Raw reads must be **paired-end Illumina FASTQ files** (gzip-compressed), organised under the working directory using the following structure:

```
<wd>/raw_reads/<sampleID>/
    <sampleID>_1.fq.gz
    <sampleID>_2.fq.gz
```

For example, for sample `sample1`:

```
/path/to/project/raw_reads/sample1/sample1_1.fq.gz
/path/to/project/raw_reads/sample1/sample1_2.fq.gz
```

### 3.2 Samples List (`samplesList.txt`)

A plain text file listing one sample ID per line, corresponding to the directory names under `raw_reads/`:

```
sample1
sample2
sample3
sample4
sample5
```

### 3.3 Samples Info (`sampleInfo.txt`)

A tab-separated file (no header) with the following columns:

| Column | Description |
|--------|-------------|
| 1 | Sample ID |
| 2 | Sex (`male` / `female`) |
| 3 | Flowcell ID |
| 4 | Lane |
| 5 | Library ID |
| 6 | Platform (e.g. `ILLUMINA`) |
| 7 | Platform unit |
| 8 | Sequencing centre |
| 9 | Affected status (`yes` = affected / `no` = control) |

Example:

```tsv
sample1	female	HK3KHDSX2	L02	lib1	ILLUMINA	ILLUMINA	ILLUMINA	yes
sample2	male	HF5L3DSX5	L04	lib1	ILLUMINA	ILLUMINA	ILLUMINA	no
sample3	female	H3FNWDSX5	L04	lib1	ILLUMINA	ILLUMINA	ILLUMINA	yes
sample4	male	HF5L3DSX5	L04	lib1	ILLUMINA	ILLUMINA	ILLUMINA	yes
sample5	male	H3FNWDSX5	L04	lib1	ILLUMINA	ILLUMINA	ILLUMINA	yes
```

---

## 4. Configuring the JSON Input File

Open `Main_inputs.json` and update the paths to match your environment. The fields that **must** be set by the user are:

```json
{
  "Main.mode": "WES",           // Analysis mode: "WES" or "WGS"
  "Main.wd": "/path/to/your/project/fastq",
  "Main.samplesList": "/path/to/your/project/fastq/raw_reads/samplesList.txt",
  "Main.samplesInfo": "/path/to/your/project/fastq/raw_reads/sampleInfo.txt",

  // --- Reference graph (pangenome) ---
  "Main.refGraph":  "./datasets/hprc-v1.1-mc-chm13.gbz",
  "Main.haplo":     "./datasets/hprc-v1.1-mc-chm13.hapl",

  // --- Linear reference paths ---
  "Main.grchPath":  "./datasets/grch38.paths.txt",
  "Main.chmPath":   "./datasets/chm13.paths.txt",
  "Main.grchDict":  "./datasets/GRCh38.dict",
  "Main.ref_hg38":  "./datasets/GRCh38_edit_chr_names.fa",
  "Main.ref_chm13": "./datasets/CHM13_edit_chr_names.fa",
  "Main.refAssembly": "./datasets/Homo_sapiens_assembly38.fasta",

  // --- Annotation and disease panel ---
  "Main.disease_gene_panel": "./datasets/disease_superpanel_example.bed",
  "Main.gnomAD":             "./datasets/gnomad.exomes.v4.0.sites.chr1-22,X,Y.vcf.gz",

  // --- Containers and helper scripts (relative paths from repo root) ---
  "Main.snv_resource":            "./containers/snvstory_resource/",
  "Main.snvstory_container":      "./containers/snvstory.sif",
  "Main.deepvariant_container":   "./containers/google_deepvariant_1.6.0.sif",
  "Main.stripy_container":        "./containers/stripy/",
  "Main.mitohpc_container":       "./containers/mitohpc.sif",
  "Main.squirls_Data":            "./containers/squirls/",
  "Main.octopanToolkit_container":"./containers/octopan-toolkit.sif",
  "Main.SnpEff":                  "./containers/snpEff",
  "Main.genderCheck":             "./containers/determine_sex_from_bam.py",
  "Main.json_to_csv_graph":       "./containers/json_to_csv_graph.R"
}
```

**Key parameter notes:**

- `Main.mode` — Set to `"WES"` for whole-exome sequencing or `"WGS"` for whole-genome sequencing. This affects variant-calling parameters and coverage expectations.
- `Main.wd` — The absolute path to the directory containing the `raw_reads/` folder.
- Relative `./datasets/` and `./containers/` paths assume you are running from the OctoPan repository root. If you run from a different directory, replace these with absolute paths.

---

## 5. Running the Pipeline

Navigate to the OctoPan repository root and execute:

```bash
java -Dconfig.file=./WDL_v3.1_fixed.conf \
     -jar ./cromwell-87.jar \
     run Main.wdl \
     -i Main_inputs.json \
     > run.log 2>&1 &
```

- Output and progress are captured in `run.log`. Monitor it with:

  ```bash
  tail -f run.log
  ```

- Cromwell will create a `cromwell-executions/` directory in the working directory containing intermediate files and per-task logs.

- On HPC clusters using SLURM, ensure `WDL_v3.1_fixed.conf` is configured to submit tasks as cluster jobs (see `cromwell_slurm.conf` for reference).

---

## 6. Workflow Steps

The pipeline proceeds through the following stages for each sample:

### Step 1 — Quality Control (`QualityControl.wdl`)
Raw FASTQ files are assessed with **FastQC** and aggregated using **MultiQC**. Reports are generated before any trimming or alignment.

### Step 2 — Pangenome Mapping (`MapToPangenome.wdl`)
Reads are aligned to the **HPRC v1.1 pangenome graph** (`hprc-v1.1-mc-chm13.gbz`) using `vg giraffe`. The output is a **GAF** (Graph Alignment Format) file, capturing alignment to the full diversity of the pangenome, including population-specific haplotypes absent from any single linear reference.

### Step 3 — BAM Surjection (`SurjectBAM.wdl`)
The pangenome GAF alignments are projected ("surjected") onto two linear reference assemblies using `vg surject`:
- **GRCh38** (hg38) — the primary clinical reference.
- **CHM13** (T2T-CHM13) — a complete telomere-to-telomere assembly with improved representation of previously unresolved regions.

Both BAMs are sorted, indexed, and made available for all downstream callers.

### Step 4 — Sex Determination (`SexDetermination.wdl`)
Chromosomal sex is inferred from coverage depth on sex chromosomes using the hg38 BAM. This result is passed to STRipy for X-linked repeat locus interpretation.

### Step 5 — SNV/Indel Variant Calling (`VariantCalling.wdl`)
**DeepVariant** (v1.6.0) calls single-nucleotide variants and small insertions/deletions independently on both hg38 and CHM13 BAMs. Variants are annotated with **SnpEff**. WES mode applies target-region constraints; WGS mode runs genome-wide.

### Step 6 — Structural Variant Calling (`Manta.wdl`)
**Manta** calls structural variants (deletions, duplications, inversions, insertions, translocations) on both hg38 and CHM13 BAMs. Running on two independent reference assemblies enables cross-validation of structural variant calls.

### Step 7 — STR Expansion Analysis (`STRipyPipeline.wdl`)
**STRipy** genotypes short-tandem repeat (STR) loci known to cause disease. Sex information from Step 4 is incorporated for correct X-linked locus interpretation. Results are produced in both JSON/CSV and HTML report formats.

### Step 8 — Splicing Variant Prioritisation (`SQUIRLS.wdl`)
**SQUIRLS** scores and prioritises splicing-relevant variants from the hg38 VCF within the supplied disease gene panel. Variants are cross-referenced against gnomAD population frequencies to filter common variants. Output is a filtered VCF and HTML report.

### Step 9 — Mutational Signature Analysis (`SNVstory.wdl`)
**SNVstory** analyses the mutational spectrum of SNVs in the hg38 VCF and decomposes them into COSMIC mutational signatures, providing insight into the aetiological processes active in the sample.

### Step 10 — Mitochondrial Extraction (`ExtractChrM.wdl`)
The mitochondrial chromosome reads are extracted from the hg38 BAM to produce a dedicated chrM BAM for mitochondrial-specific analysis.

### Step 11 — Mitochondrial Variant Analysis (`MitoHPC.wdl`)
**MitoHPC** calls mitochondrial variants, estimates heteroplasmy levels, and assigns mtDNA haplogroups. Output includes a mitochondrial VCF and a haplogroup assignment report.

### Step 12 — Tool Version Recording (`Versions.wdl`)
A manifest of all tool versions used in the run is written to a single text file for reproducibility and reporting.

---

## 7. Outputs

Cromwell writes all final outputs into the `cromwell-executions/` directory, organised by workflow name, run ID, and task name. Per-sample output paths follow the pattern:

```
cromwell-executions/Main/<run-id>/call-<TaskName>/execution/
```

The declared workflow outputs and what they contain are:

| Output variable | Description |
|-----------------|-------------|
| `FastqcReportDir` | Directory of per-sample FastQC reports |
| `MultiqcReportFile` | Aggregated MultiQC HTML summary |
| `graph` | Pangenome alignment file (GAF) |
| `bam_hg38` / `bam_hg38_idx` | Surjected BAM + index aligned to GRCh38 |
| `bam_chm13` / `bam_chm13_idx` | Surjected BAM + index aligned to CHM13 |
| `sampleSex` | Inferred chromosomal sex output |
| `vcf_file_hg38` | DeepVariant SNV/indel VCF (GRCh38) |
| `vcf_file_chm13` | DeepVariant SNV/indel VCF (CHM13) |
| `manta_vcf_file_hg38` | Manta structural variant VCF (GRCh38) |
| `manta_vcf_file_chm13` | Manta structural variant VCF (CHM13) |
| `stripy_output` | STR expansion calls (CSV/JSON) |
| `stripy_html` | STRipy HTML visual report |
| `squirls_vcf` | Splicing-prioritised variant VCF |
| `squirls_html` | SQUIRLS HTML prioritisation report |
| `snv_out` | SNVstory mutational signature output |
| `chrm_bam` | Extracted mitochondrial BAM |
| `mtdnaVCF` | Mitochondrial variant VCF |
| `haploReport` | mtDNA haplogroup assignment report |
| `toolsVersion` | Tool version manifest |

To locate a specific output after a run, inspect the Cromwell metadata or search:

```bash
find cromwell-executions/ -name "*.vcf.gz" | grep hg38
```

---

## 8. Troubleshooting

**Java not found or wrong version**
Ensure `java/openjdk-18.0.2` is installed before running. Check with `java -version`.

**R package errors in STRipy post-processing**
Install `jsonlite`, `purrr`, and `dplyr` in the R environment that Cromwell uses to execute tasks.

**Python `pandas` not found**
Install `pandas` in the Python environment available to the sex-determination script: `pip install pandas`.

**Singularity bind errors**
Ensure the project directory and `datasets/` paths are accessible to Singularity. You may need to add bind mounts to `WDL_v3.1_fixed.conf` depending on your HPC environment.

---

## Contact

For issues or questions, please open a GitHub issue
