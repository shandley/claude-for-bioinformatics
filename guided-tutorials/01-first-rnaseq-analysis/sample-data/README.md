# Sample RNA-seq Data for Tutorial

## Dataset Description
- **Type**: Paired-end RNA-seq reads
- **Read Length**: 75 bp
- **Number of Read Pairs**: 10,000
- **Total Reads**: 20,000
- **Estimated Processing Time**: 1-2 minutes for FastQC
- **File Size**: ~5-8 MB total (compressed)

## Data Characteristics
This synthetic dataset is designed for educational purposes and includes:

### Realistic Quality Patterns
- Forward reads (R1) have higher quality than reverse reads (R2)
- Quality scores decrease toward 3' end of reads (typical Illumina pattern)
- Quality scores range from ~15-40 (realistic range)
- Mix of high, medium, and some lower quality reads

### Adapter Contamination
- ~5% of reads contain adapter sequences (realistic contamination rate)
- Adapters appear at 3' end of reads with lower quality scores
- Uses common Illumina adapter sequences

### Educational Value
This dataset will produce FastQC reports that show:
- ✅ Generally good quality metrics (mostly green)
- ⚠️ Some warnings for educational discussion (yellow)
- 🔍 Realistic quality patterns for interpretation practice
- 📊 Meaningful MultiQC aggregated reports

## Expected Quality Control Results

### Per-base Quality
- Should show typical Illumina quality decline pattern
- Most bases above Q30 threshold
- Some quality drop at read ends

### Adapter Content
- Low but detectable adapter contamination (~5%)
- Good example for discussing adapter trimming

### Duplicate Levels
- Moderate duplication levels (typical for RNA-seq)
- Educational discussion about biological vs technical duplicates

### GC Content
- Approximately normal distribution
- Suitable for discussing species-specific GC content

## Usage Notes
- Generated with random seed 42 for reproducibility
- Not suitable for publication or real research
- Designed specifically for learning FastQC/MultiQC interpretation
- Processing time optimized for tutorial environment

## File Details
- `sample_R1.fastq.gz`: Forward reads (Read 1)
- `sample_R2.fastq.gz`: Reverse reads (Read 2)
- Format: Standard Illumina FASTQ with Phred+33 quality encoding
- Compression: gzip compressed for realistic file handling

Generated for Claude Code Bioinformatics Tutorial Module 1.1
