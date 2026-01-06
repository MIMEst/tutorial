# RNA-seq Analysis Pipeline  
**Gut microbiome–dependent transcriptional changes in murine colon**

This repository documents a reproducible RNA-seq analysis pipeline used to identify **gut microbiome–dependent effects on murine colon transcriptomes**, comparing **Germ-Free (GF)** and **SPF (WT)** mouse models.

## 📌 Dataset Information

- **Tissue**: Whole colon  
- **Organism**: Mouse  
- **Experimental groups**: Germ-Free vs SPF (WT)  
- **BioProject**: `PRJNA752275`  
- **Reference**:  
  Romero R, Zarzycka A, Preussner M, Fischer F *et al.*  
  *Selected commensals educate the intestinal vascular and immune system for immunocompetence.*  
  **Microbiome** (2022) 10:158  
  PMID: 36171625  

## 🧬 Pipeline Overview

1. Download raw FASTQ files from SRA  
2. Perform quality control (FastQC)  
3. Adapter trimming (Trim Galore)  
4. Alignment to reference genome (STAR)  
5. Expression quantification & normalization (RSEM)  

## ✅ Requirements

Install these tools and make sure they are in your `PATH`:

- NCBI Entrez Direct (`esearch`, `efetch`)
- SRA Toolkit (`fastq-dump` or `fasterq-dump`, `prefetch`)
- FastQC
- Trim Galore
- STAR
- RSEM

> Tip: use `conda`/`mamba` environments to manage dependencies.

---

## 1️⃣ Raw FASTQ Download

### 1.1 Retrieve SRA Run Information

SRA data are organized hierarchically as:
BioProject → BioSample → SRA Experiment → SRA Run
Download run metadata:

```bash
esearch -db sra -query PRJNA752275 | \
efetch -format runinfo > runinfo.csv
```
### 1.2 Download FASTQ Files

You may use either direct or indirect methods.
Option A: Direct download (recommended)

```bash
fastq-dump SRRXXXXXXX \
  --split-files \
  --gzip \
  --outdir fastq/
```
Option B: Indirect download via .sra (faster)
```bash
prefetch SRRXXXXXXX
fastq-dump SRRXXXXXXX/SRRXXXXXXX.sra \
  --split-files \
  --gzip \
  --outdir fastq/
```
### 2. Quality Control (FastQC)
Inspect raw FASTQ files for:
Adapter contamination
Duplication levels
Sequence length distribution
GC content


