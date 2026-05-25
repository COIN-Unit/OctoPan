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
9. [Resource Requirements and Tuning](#9-resource-requirements-and-tuning)

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

### SNPEff
```
wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
```
Unpack and store it inside the `containers/` directory.

### Squirls
```
wget https://github.com/monarch-initiative/Squirls/releases/download/v2.0.1/squirls-cli-2.0.1-distribution.zip
```
Unpack and store it inside the `containers/` directory.

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

 // --- Resources ---
  "Main.cpu":    32,     // Number of CPU cores allocated to each task
  "Main.mem_mb": 7100,   // Memory per core in MB (total RAM = cpu × mem_mb)
  "Main.time":   1200    // Maximum walltime per task in minutes
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



## 9. Resource Requirements and Tuning
 
OctoPan exposes three global resource parameters — `cpu`, `mem_mb`, and `time` — that are passed uniformly to every task in the workflow via the JSON input file. However, different tasks have very different computational demands. This section explains what each task actually needs, how to set global defaults appropriately, and when you may need to raise or lower them.
 
---
 
### 9.1 How Resources Are Controlled
 
The three parameters in `Main_inputs.json` are:
 
```json
"Main.cpu":    32,     // Number of CPU cores allocated to each task
"Main.mem_mb": 7100,   // Memory per core in MB (total RAM = cpu × mem_mb)
"Main.time":   1200    // Maximum walltime per task in minutes
```
 
> **Important:** `mem_mb` in OctoPan represents **memory per core**, not total memory. A task running with `cpu = 32` and `mem_mb = 7100` will request `32 × 7100 = 227,200 MB (~222 GB)` of RAM from the scheduler. Verify this matches your cluster's node memory limits before running.
 
These values are forwarded to Cromwell, which passes them to your HPC scheduler (e.g. SLURM) as job resource requests. If a task exceeds its allocation it will be killed by the scheduler — resulting in a failed Cromwell task.
 
---
 
### 9.2 Per-Task Resource Guide
 
The table below gives recommended settings for each workflow task for both WES (~100×) and WGS (~30×) data on a typical 150 bp paired-end Illumina library. All memory figures are **total RAM** (i.e. `cpu × mem_mb`). Adjust proportionally for higher/lower coverage or read lengths.
 
| Task | WDL module | CPU | Total RAM | Walltime (WES) | Walltime (WGS) | Notes |
|------|-----------|-----|-----------|----------------|----------------|-------|
| Quality Control | `QualityControl.wdl` | 8 | 16 GB | ~30 min | ~60 min | FastQC is per-read-file, MultiQC is lightweight. I/O bound. |
| Pangenome Mapping | `MapToPangenome.wdl` | 32 | 200–220 GB | ~3–5 h | ~8–14 h | Most memory-intensive step. `vg giraffe` loads the full GBZ graph (~60–80 GB) plus working buffers. Do not reduce CPU below 16. |
| BAM Surjection | `SurjectBAM.wdl` | 32 | 100–120 GB | ~2–3 h | ~5–8 h | `vg surject` is both CPU- and memory-intensive. Runs twice (hg38 + CHM13). |
| Sex Determination | `SexDetermination.wdl` | 4 | 8 GB | ~10 min | ~15 min | Lightweight samtools depth + Python script. |
| SNV/Indel Calling | `VariantCalling.wdl` | 32 | 120–160 GB | ~3–5 h | ~10–18 h | DeepVariant uses GPU if available; on CPU-only nodes it is the longest variant-calling step. Runs on both references. |
| Structural Variant Calling | `Manta.wdl` | 16 | 32–48 GB | ~1–2 h | ~3–6 h | Manta is highly parallel but less memory-hungry than DeepVariant. Runs on both references. |
| STR Expansion Analysis | `STRipyPipeline.wdl` | 8 | 16–24 GB | ~20–40 min | ~30–60 min | STRipy analyses a fixed set of disease loci; runtime scales little with coverage. |
| Splicing Prioritisation | `SQUIRLS.wdl` | 8 | 24–32 GB | ~30–60 min | ~60–90 min | Loads the SQUIRLS database into memory; RAM requirements are dominated by the DB size, not sample size. |
| Mutational Signatures | `SNVstory.wdl` | 8 | 16 GB | ~20–30 min | ~30–60 min | Lightweight R/Python analysis on the VCF. |
| Mitochondrial Extraction | `ExtractChrM.wdl` | 4 | 8 GB | ~5–10 min | ~5–10 min | Samtools view on a single chromosome. Very fast. |
| Mitochondrial Analysis | `MitoHPC.wdl` | 8 | 24–32 GB | ~30–60 min | ~30–60 min | Runtime is independent of nuclear genome coverage. |
| Tool Versions | `Versions.wdl` | 2 | 4 GB | ~5 min | ~5 min | Shell commands only. Negligible. |
 
---
 
### 9.3 Recommended Global Settings by Mode
 
Because `cpu`, `mem_mb`, and `time` are shared across all tasks, you must size them to the most demanding task (pangenome mapping) or accept that lighter tasks will over-request resources. The following starting points work well on most HPC nodes with ≥256 GB RAM:
 
**WES (recommended starting point)**
 
```json
"Main.cpu":    32,
"Main.mem_mb": 6500,
"Main.time":   480
```
 
Total RAM per job: `32 × 6500 = 208 GB`
Maximum walltime: 8 hours — sufficient for all WES tasks including mapping.
 
**WGS (recommended starting point)**
 
```json
"Main.cpu":    32,
"Main.mem_mb": 7100,
"Main.time":   1200
```
 
Total RAM per job: `32 × 7100 ≈ 222 GB`
Maximum walltime: 20 hours — covers the slowest WGS task (DeepVariant on both references).
 
> **Tip:** If your HPC nodes have less than 256 GB RAM, reduce `cpu` to 24 and set `mem_mb` accordingly so that `cpu × mem_mb` stays within your node's physical memory. For example: `cpu = 24`, `mem_mb = 8500` → ~200 GB total.
 
---
 
### 9.4 Memory and CPU Scaling Rules of Thumb
 
**Pangenome mapping (`vg giraffe`)** is the hardest constraint. The GBZ pangenome graph for HPRC v1.1 requires approximately 60–80 GB just to load into memory before any reads are processed. On top of that, `vg giraffe` maintains per-thread working buffers. As a rule:
 
```
minimum total RAM = graph size (~75 GB) + (cpu_threads × ~1.5 GB) + OS overhead (~10 GB)
```
 
For 32 threads: 75 + 48 + 10 = ~133 GB minimum. Allow 50–100% headroom → aim for 200+ GB.
 
**DeepVariant** requires significant RAM for model loading and per-shard processing. WGS runs are substantially longer than WES because DeepVariant processes the entire genome rather than the exome target regions. If your cluster has GPU nodes, DeepVariant can be accelerated significantly — check your Singularity/Cromwell configuration to enable GPU passthrough.
 
**Manta, STRipy, SQUIRLS, MitoHPC, SNVstory** are all substantially lighter than the mapping and calling steps. They will complete well within their allocated resources when global settings are sized for mapping/calling.
 
---
 
### 9.5 Estimating Total Pipeline Runtime
 
The pipeline runs tasks sequentially where there are data dependencies (mapping → surjection → calling) and in parallel where it can (e.g. hg38 and CHM13 calling can overlap). Approximate end-to-end wall-clock times for a single sample (assuming all tasks run immediately on a dedicated node):
 
| Mode | Optimistic | Typical | With queue wait |
|------|-----------|---------|-----------------|
| WES (100×) | ~6 h | ~10 h | ~12–24 h |
| WGS (30×) | ~14 h | ~24 h | ~36–72 h |
 
Queue wait time on shared HPC systems is environment-dependent and can dominate total turnaround time.
 
---
 
### 9.6 Setting Per-Task Resources (Advanced)
 
The current implementation uses global `cpu`, `mem_mb`, and `time` values passed uniformly to all tasks. If your Cromwell backend configuration supports per-task overrides (via the `runtime` block in individual WDL files), you can edit the relevant `.wdl` file under `wdls/` to hard-code or parameterise resources independently for each task. For example, in `MapToPangenome.wdl`:
 
```wdl
runtime {
    cpu:     cpu          # keep the global value for mapping
    memory:  mem_mb + " MB"
    time:    time
}
```
 
And in `QualityControl.wdl` you could reduce to dedicated smaller values:
 
```wdl
runtime {
    cpu:     8
    memory:  "16000 MB"
    time:    60
}
```
 
This avoids over-requesting resources for lightweight tasks on schedulers that bill by allocated (not used) resources, which can improve job priority and reduce cost on cloud HPC systems.
 
---
 
### 9.7 Monitoring Resource Usage
 
After a run completes (or fails), you can inspect actual resource consumption from the Cromwell execution logs:
 
```bash
# Check SLURM accounting for completed jobs
sacct -j <jobid> --format=JobID,JobName,MaxRSS,MaxVMSize,Elapsed,State
 
# Check peak memory from Cromwell task stderr
grep -i "memory\|oom\|killed" cromwell-executions/Main/<run-id>/call-MapToPangenomeWorkflow/execution/stderr
```
 
Use these figures to right-size `mem_mb` and `time` for future runs — particularly if jobs are being killed (OOM) or if you are wasting significant allocated-but-unused RAM.
 
---

---

## Contact

For issues or questions, please open a GitHub issue
