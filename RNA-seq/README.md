# RNA-seq Analysis Pipeline  
**Gut microbiome–dependent transcriptional changes in murine colon**

This repository presents a **reproducible bulk RNA-seq analysis pipeline** designed to investigate how the **gut microbiome shapes host transcriptional programs in the murine colon**. By comparing **Germ-Free (GF)** and **Specific Pathogen-Free (SPF; WT)** mouse models, this workflow enables systematic identification of microbiota-dependent gene expression changes at both the gene and transcript levels.

The pipeline follows widely accepted best practices in RNA-seq analysis while emphasizing **robust preprocessing, splice-aware alignment, and statistically sound expression quantification**. All steps are modular, explicitly parameterized, and suitable for reuse in other mouse or human RNA-seq studies.

---

## Dataset Information

- **Tissue**: Whole colon  
- **Organism**: Mouse (*Mus musculus*)  
- **Experimental groups**: Germ-Free vs SPF (WT)  
- **BioProject**: `PRJNA752275`  
- **Reference**:  
  Romero R, Zarzycka A, Preussner M, Fischer F *et al.*  
  *Selected commensals educate the intestinal vascular and immune system for immunocompetence.*  
  **Microbiome** (2022) 10:158  
  PMID: 36171625  

---

## Pipeline Overview and Analysis Rationale

<p align="center">
  <img width="350"
       alt="RNA-seq pipeline overview"
       src="https://github.com/user-attachments/assets/040dc06b-9d61-479b-8130-35b3adecd739" />
</p>

The analysis proceeds from raw sequencing reads to normalized expression estimates through the following conceptual stages:  
(1) acquisition of raw FASTQ files,  
(2) quality assessment of sequencing reads,  
(3) adapter and low-quality base removal,  
(4) splice-aware alignment to the reference genome, and  
(5) transcript- and gene-level expression quantification.

Each step is designed to minimize technical artifacts while preserving biologically meaningful signal, which is especially critical when comparing microbiome-driven phenotypes such as GF versus SPF intestinal tissue.

---

## Software Requirements

The following tools are required and should be available in the system `PATH`:

- NCBI Entrez Direct (`esearch`, `efetch`)
- SRA Toolkit (`prefetch`, `fastq-dump` or `fasterq-dump`)
- FastQC
- Trim Galore (Cutadapt-based)
- STAR
- RSEM

> **Tip**: Using `conda` or `mamba` environments is strongly recommended for reproducibility and version control.

---

### 1. Raw Data Retrieval and Quality Control

Raw sequencing data are retrieved from the NCBI Sequence Read Archive (SRA), which organizes data hierarchically as **BioProject → BioSample → Experiment → Run**. Run-level metadata are first obtained to identify individual sequencing accessions:

```bash
esearch -db sra -query PRJNA752275 | \
efetch -format runinfo > runinfo.csv
```

FASTQ files can then be downloaded either directly or via an intermediate .sra file. For large datasets, prefetch-based downloading is recommended:  
```bash
prefetch SRRXXXXXXX
fastq-dump SRRXXXXXXX/SRRXXXXXXX.sra \
  --split-files \
  --gzip \
  --outdir fastq/
```

Initial quality control is performed using FastQC to assess adapter contamination, base quality distributions, duplication levels, GC content, and read length profiles:  
```bash
fastqc fastq/*.fastq.gz -outdir 2.QC/
```

### 2. Adapter and Quality Trimming

Adapter sequences and low-quality bases are removed using Trim Galore, a wrapper around Cutadapt and FastQC. This step trims Illumina adapters, removes low-quality bases from read ends (default Phred < 20), and discards reads shorter than a defined minimum length, thereby improving alignment accuracy and reducing downstream noise:  
```bash
trim_galore --paired --quality 20 --length 20 \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

### 3. Splice-aware Alignment with STAR

Trimmed reads are aligned to the mouse reference genome using STAR, a splice-aware aligner optimized for RNA-seq. STAR achieves high accuracy and exceptional speed by identifying Maximal Mappable Prefixes (MMPs) through suffix-array–based searching.

<p align="center"> <img width="350" alt="STAR seed searching" src="https://github.com/user-attachments/assets/d3985f6f-e474-450a-ae2e-c7fca751185d" /> </p>

STAR iteratively maps the longest exact matches of each read (seeds) to the genome, then continues searching only the unmapped portions of the read:

<p align="center"> <img width="350" alt="STAR multiple seeds" src="https://github.com/user-attachments/assets/57ce4c8f-26b8-4ce5-b3ac-0e317f5fe5d8" /> </p>

When mismatches or indels prevent perfect matches, STAR extends alignments or applies soft clipping to low-quality or contaminating regions:

<p align="center"> <img width="350" alt="STAR extension and clipping" src="https://github.com/user-attachments/assets/94edd35d-c350-4d60-81c1-cd31ad1ff81e" /> </p>

Mapped seeds are then clustered, stitched, and scored to generate full-length spliced alignments:

<p align="center"> <img width="350" alt="STAR stitching" src="https://github.com/user-attachments/assets/6a9065db-8256-4fd9-adc2-4b71b129b7aa" /> </p>

Before alignment, the mouse reference genome and annotation are downloaded from GENCODE and indexed:  

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz
```
```bash
STAR --runThreadN 128 \
--runMode genomeGenerate \
--genomeDir ${genome_dir} \
--genomeFastaFiles ${genome.fa} \
--sjdbGTFfile ${genome.gtf}
```
Reads are then aligned and simultaneously projected into transcriptome coordinates for downstream quantification:  
```bash
STAR --runThreadN 128 \
--genomeDir ${star_index_dir} \
--readFilesIn ${sample_R1} ${sample_R2} \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outFileNamePrefix 4.align/${sample}
```

### 4. Expression Quantification and Normalization with RSEM

Expression levels are estimated using RSEM, which models read generation and assignment probabilistically to account for multi-mapping reads and transcript length bias. RSEM produces gene- and isoform-level estimates such as expected counts, TPM, and FPKM values.

<p align="center"> <img width="350" alt="RSEM workflow" src="https://github.com/user-attachments/assets/90c22704-6cb2-4393-893f-40e174aceccc" /> </p>

First, a reference is prepared using the genome and GTF annotation:

```bash
rsem-prepare-reference \
  --gtf ${gtf} \
  --STAR \
  ${genome.fa} ${rsem_index}
```

Expression is then quantified from STAR-aligned transcriptome BAM files:

```bash
rsem-calculate-expression \
  --paired-end \
  --alignments \
  -p 128 \
  ${sample_Aligned.toTranscriptome.out.bam} \
  ${rsem_index} \
  5.normalize/${sample}
```

The resulting expression matrices can be directly used for downstream differential expression analysis using tools such as DESeq2 or edgeR, enabling systematic comparison of GF and SPF transcriptional programs.


---
Good luck with your bioinformatic journey! 
