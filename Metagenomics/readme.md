# 16S rRNA-seq Analysis Pipeline  
**ASV-based gut microbiome profiling with DADA2 (paired-end Illumina)**

This repository presents a **reproducible 16S rRNA amplicon sequencing analysis pipeline** using **DADA2** in R. The workflow converts demultiplexed paired-end FASTQ files into a high-quality **ASV (amplicon sequence variant) table**, assigns taxonomy (e.g., **SILVA**), and optionally exports a **phyloseq** object for downstream microbiome analysis (alpha/beta diversity, ordination, differential abundance, etc.).

DADA2 infers **exact sequence variants** (rather than OTUs), providing higher resolution, improved reproducibility across studies, and better interpretability for strain-/variant-level signals (within the limits of 16S amplicons).

---

## Dataset Information

- **Assay**: 16S rRNA gene amplicon sequencing  
- **Reads**: Paired-end Illumina (typical 2×250 or 2×300)  
- **Input**: Demultiplexed FASTQ files (per sample)  
- **Outputs**:
  - ASV count table (`seqtab.nochim`)
  - ASV sequences (`asv_seqs.fasta`)
  - Taxonomy table (`taxa`)
  - Tracking table (`track`)
  - (Optional) `phyloseq` object (`ps.rds`)

> **Note**: This pipeline closely follows the official DADA2 tutorial:  
> https://benjjneb.github.io/dada2/tutorial.html

---

## Pipeline Overview and Analysis Rationale

<p align="center">
  <img width="520"
       alt="16S DADA2 pipeline overview"
       src="https://github.com/user-attachments/assets/your_image_id_here" />
</p>

The analysis proceeds from raw reads to ASVs through these conceptual stages:  
(1) raw read QC,  
(2) primer/adaptor removal (if needed),  
(3) quality filtering + trimming,  
(4) error model learning,  
(5) denoising and ASV inference,  
(6) paired read merging,  
(7) chimera removal,  
(8) taxonomy assignment,  
(9) downstream analysis in phyloseq.

Each step is designed to **reduce technical artifacts** (sequencing errors, chimeras, low-quality tails) while preserving true biological variation in microbial communities.

---

## Software Requirements

The following tools are required and should be available:

### R packages
- `dada2`
- `Biostrings`
- `phyloseq` (optional but recommended)
- `ggplot2` (optional)
- `readr`, `dplyr`, `tibble` (optional)

### Optional command-line tools
- `fastqc` / `multiqc` for pre-QC summaries
- `cutadapt` for primer removal (recommended if primers are still present)

> **Tip**: Use `conda/mamba` to pin versions for reproducibility.

---

## Directory Structure
├── data/  
│ ├── raw_fastq/ # input FASTQ(.gz)  
│ └── metadata/ # sample metadata (tsv/csv)  
├── results/  
│ ├── 01_qc/  
│ ├── 02_filtered/  
│ ├── 03_dada2/  
│ ├── 04_taxonomy/  
│ └── 05_phyloseq/  
├── scripts/  
│ └── 01_dada2_pipeline.R  
└── README.md  

---
## 0. Input Expectations (Important)

### FASTQ naming convention (example)
DADA2 pairs forward/reverse reads by filename. A common convention:

- `SAMPLEID_R1.fastq.gz` (forward)
- `SAMPLEID_R2.fastq.gz` (reverse)

Example:
If your files use a different pattern (e.g., `_1`/`_2`, `L001`, etc.), just adjust the `pattern=` lines in the code below.

---

## 1. (Optional) Raw QC with FastQC / MultiQC

```bash
mkdir -p results/01_qc
fastqc data/raw_fastq/*.fastq.gz -o results/01_qc
multiqc results/01_qc -o results/01_qc
```

2. DADA2 Pipeline (ASV inference)

Create and run: scripts/01_dada2_pipeline.R
```bash
## =========================================================
## 16S rRNA-seq (paired-end) DADA2 pipeline
## Based on: https://benjjneb.github.io/dada2/tutorial.html
## =========================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

# ----------------------------
# User settings
# ----------------------------
path <- "data/raw_fastq"
outdir <- "results/03_dada2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

filt_dir <- "results/02_filtered"
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)

# FASTQ patterns (edit if needed)
fnFs <- sort(list.files(path, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))

stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))

# Sample names (edit parsing if needed)
sample.names <- sub("_R1\\.fastq\\.gz$", "", basename(fnFs))

# ----------------------------
# 2.1 Inspect quality profiles
# ----------------------------
pdf(file.path(outdir, "quality_profiles.pdf"), width = 12, height = 8)
plotQualityProfile(fnFs[1:min(6, length(fnFs))])
plotQualityProfile(fnRs[1:min(6, length(fnRs))])
dev.off()

# ----------------------------
# 2.2 (Optional) Primer removal
# Recommended if primers are still present.
# Use cutadapt externally or dada2::removePrimers
# ----------------------------
# If you already removed primers, skip this section.

# ----------------------------
# 2.3 Filter and trim
# Choose truncLen based on your quality plots
# Ensure enough overlap for merging:
# overlap >= ~20 bp is a common practical minimum.
# ----------------------------
filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# Example parameters (adjust!)
# truncLen: keep high-quality region, drop low-quality tails
# maxEE: expected errors threshold
# truncQ: truncate when quality <= truncQ
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen = c(240, 200),
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

write.csv(out, file.path(outdir, "filterAndTrim_summary.csv"))

# ----------------------------
# 2.4 Learn error rates
# ----------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

pdf(file.path(outdir, "error_models.pdf"), width = 12, height = 6)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# ----------------------------
# 2.5 Dereplication
# ----------------------------
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# ----------------------------
# 2.6 Denoising (ASV inference)
# ----------------------------
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# ----------------------------
# 2.7 Merge paired reads
# ----------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# ----------------------------
# 2.8 Construct ASV table
# ----------------------------
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(outdir, "seqtab.rds"))

# ----------------------------
# 2.9 Remove chimeras
# ----------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
saveRDS(seqtab.nochim, file.path(outdir, "seqtab.nochim.rds"))

# ----------------------------
# 2.10 Track reads through pipeline
# ----------------------------
getN <- function(x) sum(getUniques(x))

track <- cbind(out,
               denoisedF = sapply(dadaFs, getN),
               denoisedR = sapply(dadaRs, getN),
               merged    = rowSums(seqtab),
               nonchim   = rowSums(seqtab.nochim))

rownames(track) <- sample.names
write.csv(track, file.path(outdir, "read_tracking.csv"))

# ----------------------------
# 2.11 Export ASV sequences (FASTA) + table (CSV)
# ----------------------------
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">ASV", seq_len(length(asv_seqs)))

fasta_file <- file.path(outdir, "asv_seqs.fasta")
writeLines(rbind(asv_headers, asv_seqs), fasta_file)

asv_tab <- t(seqtab.nochim)
colnames(asv_tab) <- sample.names
write.csv(asv_tab, file.path(outdir, "ASV_table.csv"))

message("DADA2 pipeline done.")
message("Outputs written to: ", normalizePath(outdir))
```
