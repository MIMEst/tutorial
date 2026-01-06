# RNA-seq pipeline
1. Acessing SRA and Downloading raw fastq files
* Data in the SRA is organized under a hierarchical structure
* BioProject - BioSample - SRA Experiment - SRA Run
esearch -db -sra -query PRJNA~ (bioproject) | efetch -format runinfo > runinfo.csv

2. fastq-dump
(indirect) prefetch SRR~ | fastq-dump --splite-files --gzip SRR~/SRR~.sra --out directory
(direct) fastq-dump SRR~ --split-files --gzip --out directory

3. Quality Control/Check raw fastq file
* Adapter Content, duplication level, sequence length distribution, GC% content are listed in FastQC report
fastqc illumina_100K.fastq.gz

4. Trimming adapters 
* extremely important to check which library kit is used in your sample! 
trim_galore illumina_100K.fastq.gz --illumina 
trim_galore nextera_100K.fastq.gz --nextera

5. Alignment/Mapping
