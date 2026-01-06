# RNA-seq pipeline
1. Acessing SRA and Downloading raw fastq files
* Data in the SRA is organized under a hierarchical structure
* BioProject - BioSample - SRA Experiment - SRA Run
esearch -db -sra -query PRJNA~ | efetch -format runinfo > runinfo.csv

2. fastq-dump
(indirect) prefetch SRR~ | fastq-dump --splite-files --gzip SRR~/SRR~.sra --out directory
(direct) fastq-dump SRR~ --split-files --gzip --out directory

4. Quality Control/Check raw fastq file
* Adapter Content, duplication level, sequence length distribution, GC% content are listed in FastQC report
fastqc illumina_100K.fastq.gz

4. Trimming adapters 
* Extremely important to check which library kit is used in your sample!
* automatically trimglaore search for best fit adapter! 
trim_galore illumina_100K.fastq.gz --illumina 
trim_galore nextera_100K.fastq.gz --nextera

5. Alignment/Mapping
* Be aware of where you align your processed fastq files toward!
* human/mouse genomes
- 5.1 STAR genome indexing
* genome: fasta, gtf file foramt 
(human) STAR --runThreadN number_of_threads --runMode genomeGenerate --genomeDir genome_index_output_dir --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf
- 5.2 STAR alignment 
(human) STAR --runMode alignReads --runThreadN number_of_threads --genomeDir genome_index_output_dir --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outStd Log --outFileNmaePrefix ${i} --readFilesIn read_1.gz read_2.gz --readFilesCommand zcat --quantMode TranscriptomeSAM --outSAMunmappedWithin --twopassMode Basic --outFilterMultimapNmax 1

6. Normalization
* Each tool needs individual genome indexing files
- 6.1 RSEM genome indexing 
(human) rsem-prepare-reference --gtf Homo_sapiens.GRCh38.110.gtf --STAR Homo_sapiens.GRCh38.dna.primary_assembly.fa genome_index_output_dir
- 6.2 RSEM normalizartion
(human) rsem-calculate-expression input.bam --paired-end --strandedness default --no-bam-output --alignments -p number_of_threads genome_index_output_dir


