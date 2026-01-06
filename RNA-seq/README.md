# RNA-seq Analysis Pipeline  
**Gut microbiome–dependent transcriptional changes in murine colon**

This repository documents a reproducible RNA-seq analysis pipeline used to identify **gut microbiome–dependent effects on murine colon transcriptomes**, comparing **Germ-Free (GF)** and **SPF (WT)** mouse models.

## Dataset Information

- **Tissue**: Whole colon  
- **Organism**: Mouse  
- **Experimental groups**: Germ-Free vs SPF (WT)  
- **BioProject**: `PRJNA752275`  
- **Reference**:  
  Romero R, Zarzycka A, Preussner M, Fischer F *et al.*  
  *Selected commensals educate the intestinal vascular and immune system for immunocompetence.*  
  **Microbiome** (2022) 10:158  
  PMID: 36171625  

## Pipeline Overview
<img width="300" height="700" alt="image" src="https://github.com/user-attachments/assets/040dc06b-9d61-479b-8130-35b3adecd739" />

1. Download raw FASTQ files from SRA  
2. Perform quality control (FastQC)  
3. Adapter trimming (Trim Galore)  
4. Alignment to reference genome (STAR)  
5. Expression quantification & normalization (RSEM)  

## Requirements

Install these tools and make sure they are in your `PATH`:

- NCBI Entrez Direct (`esearch`, `efetch`)
- SRA Toolkit (`fastq-dump` or `fasterq-dump`, `prefetch`)
- FastQC
- Trim Galore
- STAR
- RSEM

> Tip: use `conda`/`mamba` environments to manage dependencies.

---

## 1. Raw FASTQ Download

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
## 2. Quality Control (FastQC)
Inspect raw FASTQ files for:  
Adapter contamination  
Duplication levels  
Sequence length distribution  
GC content  
```bash
fastqc XXXXXX.fasdtq.gz -outdir 2.QC/
```
## 3. adapter trimming (Trim_galore)
* Trim_galore requires both cutadapt, fastqc for adapter trimming, Low-quality base trimming, short read removal   
* Automatically detects adapters (Illumina by default)  
* Trim_galore removes: 3′ adapter contamination, Read-through adapters (short inserts)
* Supports: Single-end, Paired-end  
```bash
trim_galore --paired --quality 20 --length 20 \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

## 4. alignment (STAR)
STAR Alignment Strategy
STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

Seed searching
Clustering, stitching, and scoring
Seed searching
For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):  
<img width="350" height="269" alt="image" src="https://github.com/user-attachments/assets/d3985f6f-e474-450a-ae2e-c7fca751185d" />

The different parts of the read that are mapped separately are called ‘seeds’. So the first MMP that is mapped to the genome is called seed1.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be seed2.  
<img width="350" height="204" alt="image" src="https://github.com/user-attachments/assets/57ce4c8f-26b8-4ce5-b3ac-0e317f5fe5d8" />

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

If STAR does not find an exact matching sequence for each part of the read due to mismatches or indels, the previous MMPs will be extended.  
<img width="350" height="460" alt="image" src="https://github.com/user-attachments/assets/94edd35d-c350-4d60-81c1-cd31ad1ff81e" />

If extension does not give a good alignment, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.  
<img width="350" height="213" alt="image" src="https://github.com/user-attachments/assets/4a74e78b-117b-4820-89f1-0b725628bd89" />

Clustering, stitching, and scoring
The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).  
<img width="350" height="288" alt="image" src="https://github.com/user-attachments/assets/6a9065db-8256-4fd9-adc2-4b71b129b7aa" />

### Download mouse genome and gtf file 
GENCODE only provides human and mouse annotations. The genome and annotation files can be found from GENCODE website.
Download and decompress the mouse genome and GTF files:  
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz
```

### 4.1 Creat a genome index
building mouse/human genome index  
(!please be aware of your dataset!)
```bash
STAR --runThreadN 128 \
--runMode genomeGenerate \
--genomeDir ${genome_dir} \
--genomeFastaFiles ${genome.fa} \
--sjdbGTFfile ${genome.gtf} \
--outdir ${output_dir}
```

### 4.2 Map reads to the genome
After generating genome index, align paired reads into mouse/human genome
```bash
STAR --runMode alignReads \
--runThreadN 128 \
--genomeDir ${star_index_dir} \
--outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout \
--outStd Log --outFileNamePrefix 4.align/${sample} \
--readFilesIn ${sample_R1} ${sample_R2} \
--readFilesCommand zcat \
--quantMode TranscriptomeSAM --limitBAMsortRAM 10000000000
```

## 5. Normalization (RSEM)
The rsem-calculate-expression script handles both the alignment of reads against reference transcript sequences and the calculation of relative abundances. By default, RSEM uses the Bowtie alignment program to align reads, with parameters specifically chosen for RNA-Seq quantification. Alternatively, users may manually run a different alignment program and provide alignments in SAM format to rsem-calculate-expression  
<img width="300" height="700" alt="image" src="https://github.com/user-attachments/assets/90c22704-6cb2-4393-893f-40e174aceccc" />  

RSEM is also two-step process proceeding genome indexing and normalization

### 5.1 Create genome index
Preparing Reference Sequences
```bash
rsem-prepare-reference --gtf ${gtf} --STAR ${genome.fa} ${output_dir}
```

### 5.2 Estimate gene and isoform expression
Calculating Expression Values
```bash
rsem-calculate-expression --paired-end --no-bam-output --alignments -p 128 ${sample_Aligned.toTranscriptome.out.bam} \
${rsem_genome_index} 5.normalize/${sample}
```




