# Bioinformatics Computational Resources Guide

## Table of Contents
1. [Memory Requirements by Tool](#memory-requirements-by-tool)
2. [Runtime Estimates](#runtime-estimates)
3. [Storage Planning](#storage-planning)
4. [Parallelization Strategies](#parallelization-strategies)
5. [HPC Job Submission](#hpc-job-submission)
6. [Performance Optimization](#performance-optimization)
7. [Resource Monitoring](#resource-monitoring)

---

## Memory Requirements by Tool

### Sequence Alignment Tools

#### BWA-MEM
```bash
# Memory usage: ~3GB + reference genome size
# Human genome (GRCh38): ~6-8GB total
# Mouse genome: ~4-5GB total
# Bacterial genome: ~3-4GB total

# Conservative estimates:
# Small genome (<100MB): 4GB
# Medium genome (500MB-3GB): 8GB  
# Large genome (>3GB): 16GB

# Example job submission:
#SBATCH --mem=8G  # For human genome
bwa mem -t 8 reference.fa reads_R1.fq reads_R2.fq
```

#### STAR Alignment
```bash
# Memory usage: Genome index size + 1-2GB overhead
# Human genome index: ~28-32GB
# Mouse genome index: ~25-28GB
# Drosophila genome index: ~4-6GB

# Memory by species:
# Human: 32-64GB (depending on annotation complexity)
# Mouse: 28-32GB
# Zebrafish: 16-20GB
# Drosophila: 8-12GB
# C. elegans: 4-8GB
# E. coli: 2-4GB

# Example job submission:
#SBATCH --mem=32G  # For human genome
STAR --runThreadN 8 --genomeDir star_index --readFilesIn R1.fq R2.fq
```

#### HISAT2
```bash
# Memory usage: More efficient than STAR
# Human genome: 8-12GB
# Mouse genome: 6-10GB
# Smaller genomes: 4-8GB

# Example:
#SBATCH --mem=12G
hisat2 -x genome_index -1 R1.fq -2 R2.fq -p 8
```

#### minimap2
```bash
# Memory usage: Very efficient
# Human genome: 4-8GB
# Long-read alignment: +2GB per Gb of reads

# Example:
#SBATCH --mem=8G
minimap2 -t 8 -ax map-ont reference.fa longreads.fq
```

### Variant Calling Tools

#### GATK HaplotypeCaller
```bash
# Memory per sample:
# WGS: 8-16GB per sample
# WES: 4-8GB per sample
# Targeted panels: 2-4GB per sample

# Joint genotyping scales with sample number:
# 10 samples: 16-32GB
# 100 samples: 64-128GB
# 1000+ samples: Use GenomicsDB (32-64GB)

# Example:
#SBATCH --mem=16G  # Single sample WGS
gatk --java-options "-Xmx12g" HaplotypeCaller -R ref.fa -I input.bam
```

#### FreeBayes
```bash
# Memory usage: Generally lower than GATK
# WGS: 4-8GB per sample
# Scales linearly with sample number

# Example:
#SBATCH --mem=8G
freebayes -f reference.fa input.bam > variants.vcf
```

#### DeepVariant
```bash
# Memory usage: GPU-dependent
# CPU mode: 16-32GB
# GPU mode: 8-16GB + GPU memory

# Example:
#SBATCH --mem=24G --gres=gpu:1
run_deepvariant --model_type=WGS --ref=ref.fa --reads=input.bam
```

### Assembly Tools

#### SPAdes
```bash
# Memory scaling with data size:
# Bacterial genome (5MB): 8-16GB
# Small eukaryote (100MB): 32-64GB  
# Large genome (1GB+): 128-512GB

# Rule of thumb: 50-100x the genome size
# For bacterial genomes: 16-32GB usually sufficient
# For mammalian genomes: 256-512GB required

# Example bacterial assembly:
#SBATCH --mem=32G
spades.py -1 R1.fq -2 R2.fq -o output_dir --careful
```

#### Canu
```bash
# Memory usage for long-read assembly:
# Bacterial genome: 16-32GB
# Mammalian genome: 64-256GB
# Large plant genomes: 256-512GB

# Example:
#SBATCH --mem=64G
canu -p assembly -d output genomeSize=5m -pacbio-raw reads.fq
```

#### Flye
```bash
# Memory usage: More efficient than Canu
# Bacterial: 8-16GB
# Mammalian: 32-128GB

# Example:
#SBATCH --mem=32G
flye --pacbio-raw reads.fq --out-dir output --genome-size 5m
```

### RNA-seq Analysis Tools

#### Salmon
```bash
# Memory usage: Very efficient
# Human transcriptome: 4-8GB
# Most organisms: 2-4GB

# Example:
#SBATCH --mem=8G
salmon quant -i index -l A -1 R1.fq -2 R2.fq -p 8
```

#### DESeq2 (R)
```bash
# Memory depends on sample number and gene count:
# 10 samples, 20k genes: 4-8GB
# 100 samples, 20k genes: 16-32GB
# 1000 samples, 20k genes: 64-128GB

# Example R session:
# Request 32GB for medium datasets
options(java.parameters = "-Xmx30g")
```

### Single-Cell Analysis

#### Cell Ranger
```bash
# Memory usage by cell number:
# 1K-10K cells: 16-32GB
# 10K-100K cells: 64-128GB
# 100K+ cells: 256-512GB

# Example:
#SBATCH --mem=64G
cellranger count --id=sample --transcriptome=refdata --fastqs=fastq_path
```

#### scanpy (Python)
```bash
# Memory scales with cell and gene number:
# 10K cells: 8-16GB
# 100K cells: 32-64GB
# 1M cells: 128-256GB

# Use sparse matrices for efficiency
import scanpy as sc
sc.settings.max_memory = 64  # GB
```

---

## Runtime Estimates

### Alignment Runtimes

#### Human WGS (30x coverage, ~100GB FASTQ)
```bash
# BWA-MEM (8 cores): 6-12 hours
# STAR (8 cores): 2-4 hours
# minimap2 (8 cores): 1-2 hours
# HISAT2 (8 cores): 3-6 hours

# Factors affecting runtime:
# - Read length (longer = slower)
# - Read quality (poor quality = slower)
# - Reference complexity (more repeats = slower)
```

#### Human WES (100x coverage, ~20GB FASTQ)
```bash
# BWA-MEM (8 cores): 2-4 hours
# STAR (8 cores): 1-2 hours
# HISAT2 (8 cores): 1-3 hours
```

#### RNA-seq (50M paired-end reads)
```bash
# STAR (8 cores): 30-60 minutes
# HISAT2 (8 cores): 45-90 minutes
# Salmon (8 cores): 5-15 minutes
```

### Variant Calling Runtimes

#### GATK HaplotypeCaller
```bash
# Human WGS (30x): 12-24 hours (single core)
# Human WGS (30x): 4-8 hours (8 cores with scatter-gather)
# Human WES (100x): 2-6 hours (single core)

# Joint genotyping:
# 10 samples: 2-4 hours
# 100 samples: 8-16 hours
# 1000 samples: Use GenomicsDB workflow
```

#### FreeBayes
```bash
# Generally 2-3x faster than GATK HaplotypeCaller
# Human WGS: 4-8 hours (8 cores)
# Human WES: 1-3 hours (8 cores)
```

### Assembly Runtimes

#### Bacterial Genome Assembly
```bash
# SPAdes (30x coverage): 1-4 hours
# Canu (PacBio, 50x): 4-8 hours
# Flye (PacBio, 50x): 2-4 hours
```

#### Mammalian Genome Assembly
```bash
# SPAdes (not recommended for large genomes)
# Canu (PacBio): 2-7 days
# Flye (PacBio): 1-3 days
```

---

## Storage Planning

### Data Size Estimates

#### Raw Sequencing Data
```bash
# Human WGS (30x coverage):
# FASTQ: 90-120GB per sample
# Compressed: 30-40GB per sample

# Human WES (100x coverage):
# FASTQ: 20-30GB per sample
# Compressed: 8-12GB per sample

# RNA-seq (50M reads):
# FASTQ: 15-25GB per sample
# Compressed: 5-8GB per sample

# Single-cell RNA-seq (10K cells):
# FASTQ: 5-15GB per sample
# Compressed: 2-5GB per sample
```

#### Processed Data
```bash
# BAM files (aligned reads):
# WGS: 30-50GB per sample
# WES: 5-10GB per sample
# RNA-seq: 3-8GB per sample

# VCF files:
# Single sample WGS: 1-5GB
# Cohort (100 samples): 10-50GB
# Large population: Use BCF format

# Assembly outputs:
# Bacterial genome: 5-10MB
# Mammalian genome: 2-4GB
```

#### Intermediate Files
```bash
# Plan for 2-3x raw data size in temporary storage
# GATK workflows can generate large temporary files
# Assembly can require 10-20x genome size in temp space

# Example for human WGS:
# Raw FASTQ: 100GB
# Intermediate files: 200-300GB
# Final outputs: 50GB
# Total peak usage: 400GB per sample
```

### Storage Architecture

#### Recommended Layout
```bash
/project/
├── data/
│   ├── raw/           # Raw sequencing data
│   ├── processed/     # Quality-controlled data
│   └── reference/     # Reference genomes and indices
├── analysis/
│   ├── alignment/     # BAM files
│   ├── variants/      # VCF files
│   └── expression/    # Gene counts, etc.
├── results/           # Final outputs
├── scripts/           # Analysis scripts
└── scratch/           # Temporary files (auto-delete)
```

#### Backup Strategy
```bash
# Critical data (irreplaceable):
# - Raw sequencing data
# - Final analysis results
# - Analysis scripts and metadata

# Can be regenerated (consider not backing up):
# - Intermediate alignment files
# - Temporary analysis files
# - Large intermediate assemblies
```

---

## Parallelization Strategies

### Thread-Based Parallelization

#### Optimal Thread Counts
```bash
# CPU-bound tools (alignment, assembly):
# Use number of physical cores
# Avoid hyperthreading for CPU-intensive tasks

# I/O-bound tools (file processing):
# Can use 1.5-2x physical cores

# Memory-bound tools:
# Fewer threads to avoid memory pressure

# Example scaling:
# 8-core system: use 6-8 threads for CPU-bound
# 16-core system: use 12-16 threads for CPU-bound
# 32-core system: use 24-32 threads (watch memory)
```

#### Tool-Specific Threading
```bash
# BWA-MEM: scales well to 8-16 threads
bwa mem -t 8 ref.fa reads.fq

# STAR: scales well to 8-32 threads
STAR --runThreadN 16

# GATK: most tools are single-threaded
# Use scatter-gather for parallelization
gatk --java-options "-XX:ParallelGCThreads=4"

# samtools: most operations scale linearly
samtools sort -@ 8 input.bam
```

### Data Parallelization

#### Chromosome-Based Splitting
```bash
# Process each chromosome separately
for chr in {1..22} X Y; do
    gatk HaplotypeCaller -L chr${chr} -I input.bam &
done
wait

# Combine results
gatk GatherVcfs -I chr*.vcf -O combined.vcf
```

#### Region-Based Splitting
```bash
# Split genome into regions
gatk SplitIntervals -R ref.fa --scatter-count 32 -O intervals/

# Process each interval
for interval in intervals/*.interval_list; do
    gatk HaplotypeCaller -L $interval -I input.bam &
done
```

#### Sample-Based Parallelization
```bash
# Process multiple samples in parallel
for sample in sample*.bam; do
    gatk HaplotypeCaller -I $sample &
    # Limit concurrent jobs
    (($(jobs -r | wc -l) >= 4)) && wait
done
```

---

## HPC Job Submission

### SLURM Job Scripts

#### Single-Sample Analysis
```bash
#!/bin/bash
#SBATCH --job-name=wgs_analysis
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=analysis_%j.log
#SBATCH --error=analysis_%j.err

# Load modules
module load bwa samtools gatk

# Set temporary directory
export TMPDIR=/scratch/$USER/$SLURM_JOB_ID
mkdir -p $TMPDIR

# Analysis commands
bwa mem -t 8 ref.fa R1.fq R2.fq | samtools sort -@ 4 -o aligned.bam
gatk --java-options "-Xmx28g" HaplotypeCaller -R ref.fa -I aligned.bam

# Cleanup
rm -rf $TMPDIR
```

#### Array Job for Multiple Samples
```bash
#!/bin/bash
#SBATCH --job-name=multi_sample
#SBATCH --array=1-100
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Get sample name from array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

# Process sample
bwa mem -t 4 ref.fa ${SAMPLE}_R1.fq ${SAMPLE}_R2.fq | \
    samtools sort -@ 2 -o ${SAMPLE}.bam
```

#### Dependency Jobs
```bash
# Submit alignment job
ALIGN_JOB=$(sbatch align.sh | awk '{print $4}')

# Submit variant calling dependent on alignment
sbatch --dependency=afterok:$ALIGN_JOB variant_calling.sh
```

### PBS/Torque Job Scripts

#### Basic PBS Script
```bash
#!/bin/bash
#PBS -N analysis_job
#PBS -l walltime=24:00:00
#PBS -l mem=32gb
#PBS -l ncpus=8
#PBS -o analysis.log
#PBS -e analysis.err

cd $PBS_O_WORKDIR
module load biotools

# Analysis commands here
```

### Resource Estimation Guidelines

#### Conservative Estimates
```bash
# Memory: Request 1.5-2x expected usage
# Time: Request 1.5-2x expected runtime
# CPUs: Start with 4-8, scale based on tool efficiency

# Example for human WGS alignment:
# Expected: 8GB memory, 6 hours, 8 CPUs
# Request: 16GB memory, 12 hours, 8 CPUs
```

#### Monitoring and Adjustment
```bash
# Check resource usage after jobs complete
sacct -j JOBID --format=JobID,MaxRSS,Elapsed,CPUTime
seff JOBID  # Efficiency report

# Adjust future jobs based on actual usage
```

---

## Performance Optimization

### I/O Optimization

#### Use Local Storage
```bash
# Copy data to local scratch
cp /network/data/* /scratch/$USER/
# Run analysis on local storage
# Copy results back
cp results/* /network/output/
```

#### Streaming Operations
```bash
# Avoid intermediate files
bwa mem ref.fa reads.fq | samtools sort -o output.bam

# Use pipes for multi-step operations
samtools view input.bam | \
awk '{if($3 != "*") print}' | \
samtools view -bS - > filtered.bam
```

#### Compression Strategies
```bash
# Use appropriate compression
# FASTQ: gzip or better yet, use --readFilesCommand zcat
# BAM: automatically compressed
# VCF: use bgzip for indexing compatibility

# Parallel compression
pigz -p 8 large_file.fastq  # Faster than gzip
```

### Memory Optimization

#### Efficient Data Structures
```bash
# Use sparse matrices for single-cell data
# Use streaming readers for large files
# Implement garbage collection in long scripts

# Python example:
import gc
# After large data processing
del large_object
gc.collect()
```

#### Memory Mapping
```bash
# For tools that support it
samtools view -M input.bam  # Memory mapping

# For custom scripts, use memory-mapped files
# Python: numpy.memmap
# R: bigmemory package
```

---

## Resource Monitoring

### Real-Time Monitoring

#### System Resources
```bash
# Memory usage
free -h
cat /proc/meminfo | grep Available

# CPU usage
top -u $USER
htop

# I/O usage
iotop
iostat 1

# GPU usage (if applicable)
nvidia-smi -l 1
```

#### Process-Specific Monitoring
```bash
# Monitor specific process
pid=$(pgrep tool_name)
while kill -0 $pid 2>/dev/null; do
    ps -p $pid -o pid,pcpu,pmem,etime,cmd
    sleep 60
done

# Memory usage over time
while true; do
    echo "$(date): $(ps -p $pid -o rss= | tail -n1) KB"
    sleep 300
done > memory_usage.log
```

### Job Accounting

#### SLURM Accounting
```bash
# Check job efficiency
seff JOBID

# Detailed accounting
sacct -j JOBID --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode

# Resource usage over time
sacct -j JOBID --format=JobID,MaxRSS,MaxVMSize,AveCPU,Elapsed
```

#### Custom Monitoring Scripts
```bash
#!/bin/bash
# Monitor resource usage during analysis

LOGFILE="resource_usage_$(date +%Y%m%d_%H%M%S).log"

while true; do
    echo "$(date)" >> $LOGFILE
    echo "Memory:" >> $LOGFILE
    free -h >> $LOGFILE
    echo "Disk:" >> $LOGFILE
    df -h >> $LOGFILE
    echo "Processes:" >> $LOGFILE
    ps aux | grep $USER | head -20 >> $LOGFILE
    echo "---" >> $LOGFILE
    sleep 300
done
```

---

## Cost Optimization

### Cloud Computing Considerations

#### Instance Selection
```bash
# CPU-optimized instances for alignment
# Memory-optimized instances for assembly
# GPU instances for deep learning tools

# Spot instances for non-critical jobs
# Reserved instances for long-running analyses
```

#### Data Transfer Optimization
```bash
# Minimize data movement
# Use compression for network transfers
# Consider data locality in cloud regions

# Example: compress before upload
tar czf data.tar.gz data/
aws s3 cp data.tar.gz s3://bucket/
```

### On-Premise Optimization

#### Queue Management
```bash
# Use appropriate job queues
# Short jobs: quick queue (< 1 hour)
# Long jobs: long queue (> 24 hours)
# High memory jobs: bigmem queue

# Example SLURM partitions
#SBATCH --partition=quick    # For short jobs
#SBATCH --partition=long     # For long jobs
#SBATCH --partition=bigmem   # For memory-intensive jobs
```

#### Resource Sharing
```bash
# Use shared reference databases
# Implement workflow caching
# Share intermediate results when possible

# Example: shared reference location
REF_DIR="/shared/references/GRCh38"
export BWA_INDEX="$REF_DIR/bwa_index/GRCh38"
export STAR_INDEX="$REF_DIR/star_index"
```

---

*This computational resources guide provides practical guidelines for planning and optimizing bioinformatics analyses across different computing environments. Resource requirements may vary based on specific implementations and data characteristics.*