# Bioinformatics Troubleshooting Guide

## Table of Contents
1. [Memory and Resource Errors](#memory-and-resource-errors)
2. [File Format Issues](#file-format-issues)
3. [Tool-Specific Error Messages](#tool-specific-error-messages)
4. [Environment and Dependencies](#environment-and-dependencies)
5. [Data Quality Issues](#data-quality-issues)
6. [Performance Optimization](#performance-optimization)
7. [Common Workflow Failures](#common-workflow-failures)

---

## Memory and Resource Errors

### Out of Memory (OOM) Errors

#### Symptoms
- "Killed" messages in logs
- Exit code 137
- Process terminates unexpectedly
- System becomes unresponsive

#### Common Tools Affected
- STAR alignment (genome loading)
- BWA-MEM with large genomes
- GATK joint genotyping
- Assembly tools (SPAdes, Canu)

#### Solutions by Tool

**STAR Alignment:**
```bash
# Problem: Insufficient memory for genome loading
# Solution: Use shared memory or reduce genome
STAR --genomeLoad LoadAndKeep --limitGenomeGenerateRAM 31000000000
# Or use smaller genome index
STAR --genomeSAindexNbases 13  # For smaller genomes
```

**BWA-MEM:**
```bash
# Problem: Large genome + many reads
# Solution: Process in smaller batches
split -l 4000000 large_file.fastq split_  # Split FASTQ
for file in split_*; do
    bwa mem ref.fa $file > ${file}.sam
done
```

**GATK HaplotypeCaller:**
```bash
# Problem: Large cohort joint genotyping
# Solution: Use GenomicsDB for scalability
gatk GenomicsDBImport --batch-size 50
gatk GenotypeGVCFs --genomicsdb-workspace
```

#### General Memory Optimization
- **Monitor usage**: `htop` or `top` during analysis
- **Estimate requirements**: Rule of thumb: 3x data size for alignment
- **Use streaming**: Pipe commands to avoid intermediate files
- **Split data**: Process chromosomes or regions separately
- **Increase swap**: Temporary solution for small overruns

### Disk Space Issues

#### Symptoms
- "No space left on device"
- Write errors during analysis
- Incomplete output files

#### Solutions
```bash
# Check disk usage
df -h
du -sh * | sort -hr

# Clean temporary files
find . -name "*.tmp" -delete
find . -name "core.*" -delete

# Compress large files
gzip *.fastq *.sam *.vcf
pigz -p 8 large_file.fastq  # Parallel compression

# Use streaming to avoid intermediate files
bwa mem ref.fa R1.fq R2.fq | samtools sort -o output.bam
```

---

## File Format Issues

### Corrupt or Truncated Files

#### Detection
```bash
# Check FASTQ integrity
seqkit stats file.fastq.gz  # Should show even number of records
zcat file.fastq.gz | tail -n 4  # Should end with quality line

# Check BAM integrity
samtools quickcheck file.bam
samtools view file.bam | tail -n 1  # Should show complete record

# Check VCF integrity
bcftools view file.vcf.gz | tail -n 5
```

#### Common Fixes
```bash
# Fix truncated FASTQ
seqkit seq --validate-seq file.fastq > cleaned.fastq

# Repair BAM files
samtools view -h broken.bam | samtools view -bS - > fixed.bam

# Fix VCF headers
bcftools view -h broken.vcf > header.txt
# Edit header manually, then:
bcftools reheader -h header.txt broken.vcf > fixed.vcf
```

### Encoding Issues

#### Phred Quality Score Problems
```bash
# Detect encoding
seqkit stats file.fastq  # Shows quality encoding
# Phred+33: range 33-126 (Sanger)
# Phred+64: range 64-126 (Illumina 1.3+)

# Convert Phred+64 to Phred+33
seqtk seq -Q64 -V input.fq > output.fq
```

#### Line Ending Issues
```bash
# Fix Windows line endings
dos2unix file.txt
sed -i 's/\r$//' file.txt

# Fix Mac line endings
sed -i 's/\r/\n/g' file.txt
```

### Chromosome Naming Inconsistencies

#### Problem: Reference vs data chromosome naming
```bash
# Check chromosome names
samtools view file.bam | cut -f3 | sort | uniq
bcftools view file.vcf | grep -v "^#" | cut -f1 | sort | uniq

# Common mismatches:
# - "chr1" vs "1"
# - "chrM" vs "MT"
# - Different scaffold naming

# Fix chromosome names in BAM
samtools view -H file.bam | sed 's/SN:chr/SN:/' | samtools reheader - file.bam > fixed.bam

# Fix chromosome names in VCF
bcftools annotate --rename-chrs chr_name_map.txt input.vcf.gz > output.vcf.gz
```

---

## Tool-Specific Error Messages

### GATK Errors

#### "A USER ERROR has occurred"
**Common Causes & Solutions:**

```bash
# Error: Input files have different sequence dictionaries
# Solution: Ensure all files use same reference
gatk CreateSequenceDictionary -R reference.fa

# Error: Read group information is missing
# Solution: Add read groups during alignment
bwa mem -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA" ref.fa reads.fq > output.sam

# Error: Unmapped reads in file
# Solution: Filter out unmapped reads
samtools view -F 4 -b input.bam > mapped_only.bam
```

#### GATK Memory Errors
```bash
# Error: GC overhead limit exceeded
# Solution: Increase heap size
gatk --java-options "-Xmx8g" HaplotypeCaller

# Error: OutOfMemoryError
# Solution: Use streaming or split regions
gatk HaplotypeCaller -L chr1:1-50000000  # Process in chunks
```

### BWA Errors

#### "Inconsistent mate pair information"
```bash
# Solution: Fix paired-end reads
samtools fixmate input.bam fixed.bam
samtools sort fixed.bam > sorted.bam
```

#### "Different number of sequences"
```bash
# Error: R1 and R2 files have different read counts
# Solution: Repair paired files
repair.sh in1=R1.fq in2=R2.fq out1=fixed_R1.fq out2=fixed_R2.fq
```

### STAR Errors

#### "EXITING because of FATAL ERROR in reads input"
```bash
# Error: Inconsistent read lengths
# Solution: Trim to consistent length
cutadapt --length 100 input.fastq > trimmed.fastq

# Error: Mate pairs are not consistent
# Solution: Sort paired files
sortByName.sh in1=R1.fq in2=R2.fq out1=sorted_R1.fq out2=sorted_R2.fq
```

#### "EXITING because of FATAL ERROR: could not create output file"
```bash
# Error: Permission or space issues
# Solution: Check directory permissions and space
chmod 755 output_directory/
df -h output_directory/
```

### samtools/bcftools Errors

#### "[E::hts_open_format] fail to open file"
```bash
# Check file exists and permissions
ls -la file.bam
file file.bam  # Should show "data" not "text"

# Check if file is corrupted
samtools quickcheck file.bam
```

#### "invalid BAM binary header"
```bash
# Rebuild BAM header
samtools view input.bam | samtools view -bT reference.fa - > fixed.bam
```

---

## Environment and Dependencies

### Conda/Mamba Environment Issues

#### Package Conflicts
```bash
# Check for conflicts
conda list | grep -E "(gatk|bwa|star)"

# Create clean environment
conda create -n bioinf python=3.9
conda activate bioinf
conda install -c bioconda gatk4 bwa star samtools

# Or use mamba for faster solving
mamba install -c bioconda tool_name
```

#### Library Version Conflicts
```bash
# Error: GLIBC version conflicts
# Solution: Use bioconda packages or containers
conda install -c bioconda glibc

# Check library dependencies
ldd $(which tool_name)
```

### PATH and Module Issues

#### Command Not Found
```bash
# Check if tool is installed
which tool_name
echo $PATH

# For HPC systems with modules
module avail bio
module load biotools
module list

# Add to PATH temporarily
export PATH=/path/to/tool:$PATH
```

#### Java Version Issues
```bash
# GATK requires Java 8+
java -version
export JAVA_HOME=/path/to/java8

# Or use conda java
conda install openjdk=8
```

---

## Data Quality Issues

### Poor Sequence Quality

#### Low Quality Scores
```bash
# Assess quality
fastqc input.fastq.gz
multiqc .

# Trim poor quality
trimmomatic PE input_R1.fq input_R2.fq \
    output_R1.fq unpaired_R1.fq \
    output_R2.fq unpaired_R2.fq \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Or use fastp
fastp -i input.fq -o output.fq -q 20 -l 30
```

### Adapter Contamination

#### Detection and Removal
```bash
# Detect adapters
fastqc input.fastq.gz  # Check "Adapter Content" section

# Remove adapters
cutadapt -a AGATCGGAAGAG -o trimmed.fq input.fq
# For paired-end
cutadapt -a ADAPT1 -A ADAPT2 -o out1.fq -p out2.fq in1.fq in2.fq
```

### Duplication Issues

#### High PCR Duplication
```bash
# Mark duplicates
gatk MarkDuplicates -I input.bam -O marked.bam -M metrics.txt

# Check duplication rate
samtools flagstat marked.bam
# Look for "duplicates" percentage

# If >30%, consider:
# - Library complexity issues
# - Over-amplification during PCR
# - Optical duplicates (check cluster distance)
```

---

## Performance Optimization

### Slow Analysis Performance

#### General Optimization Strategies
```bash
# Use multiple threads
bwa mem -t 8 ref.fa reads.fq > output.sam
samtools sort -@ 8 input.bam > sorted.bam

# Use faster tools when appropriate
minimap2 -t 8 ref.fa reads.fq > output.sam  # Faster than BWA for long reads

# Pipeline optimization
bwa mem ref.fa reads.fq | samtools sort -@ 4 -o output.bam
```

#### I/O Optimization
```bash
# Use local storage for temporary files
export TMPDIR=/local/scratch
samtools sort -T /local/scratch/temp input.bam > output.bam

# Compress intermediate files
bgzip large_file.vcf
tabix -p vcf large_file.vcf.gz
```

### Network File System Issues

#### Slow NFS Performance
```bash
# Copy to local storage first
cp /nfs/data/* /local/scratch/
# Run analysis on local data
# Copy results back
cp results/* /nfs/output/
```

---

## Common Workflow Failures

### Pipeline Crashes

#### Identify Failure Points
```bash
# Check exit codes
echo $?  # Non-zero indicates error

# Check log files
tail -n 50 analysis.log
grep -i error analysis.log
grep -i warning analysis.log

# Check intermediate files
ls -la *.bam *.vcf
samtools quickcheck *.bam
```

#### Recovery Strategies
```bash
# Resume from checkpoint
# Most pipelines support --resume or --continue flags

# Restart from specific step
nextflow run pipeline.nf --resume
snakemake --rerun-incomplete

# Manual recovery
# Identify last successful step
# Remove failed outputs
# Restart from that point
```

### Resource Allocation Failures

#### HPC Job Failures
```bash
# Check job status
squeue -u $USER
sacct -j JOBID

# Common issues:
# - Insufficient memory requested
# - Time limit exceeded
# - Node failure

# Adjust resource requests
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
```

---

## Quick Diagnostic Commands

### System Resources
```bash
# Memory usage
free -h
cat /proc/meminfo

# CPU usage
top
htop

# Disk usage
df -h
du -sh *

# Process monitoring
ps aux | grep tool_name
pgrep -f tool_name
```

### File Validation
```bash
# Quick file checks
file *.bam *.vcf.gz
wc -l *.fastq
md5sum file.txt  # Check integrity

# Format validation
samtools quickcheck *.bam
bcftools view -h file.vcf.gz | head
seqkit stats *.fastq.gz
```

### Network and Connectivity
```bash
# Check network access
ping google.com
wget -q --spider http://example.com

# Check database access
wget -q --spider ftp://ftp.ncbi.nlm.nih.gov/
```

---

## Emergency Recovery Procedures

### Data Recovery
```bash
# Recover from backups
rsync -av backup_location/ ./

# Partial file recovery
dd if=corrupted_file of=recovered_file bs=512 conv=noerror,sync

# Check for auto-saves or temporary files
find . -name "*.tmp" -o -name "*~" -o -name "*.bak"
```

### Process Recovery
```bash
# Kill runaway processes
pkill -f tool_name
killall tool_name

# Clean up shared memory (STAR)
ipcs -m | grep $USER | awk '{print $2}' | xargs -I{} ipcrm -m {}

# Reset terminal
reset
stty sane
```

---

*This troubleshooting guide covers the most common issues encountered in bioinformatics workflows. For tool-specific errors not covered here, consult the tool's documentation or check the project's GitHub issues page.*