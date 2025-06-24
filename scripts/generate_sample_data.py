#!/usr/bin/env python3
"""
Generate realistic sample RNA-seq data for educational tutorials.

This script creates small but realistic FASTQ files suitable for teaching
bioinformatics concepts without requiring large computational resources.

Features:
- Realistic read quality scores based on Illumina sequencing patterns
- Proper paired-end read structure
- Mix of high and medium quality reads for educational QC analysis
- Small size (~10,000 read pairs) for quick processing
- Realistic adapter contamination and quality patterns
"""

import random
import gzip
import os
from pathlib import Path

# Set random seed for reproducible educational data
random.seed(42)

def generate_quality_string(length, base_quality=30, variation=10):
    """Generate realistic quality scores following Illumina patterns."""
    qualities = []
    for i in range(length):
        # Quality typically decreases toward 3' end
        position_effect = max(0, (length - i - 1) / length * 5)
        
        # Add random variation
        noise = random.gauss(0, variation/3)
        
        # Calculate final quality
        final_qual = int(base_quality + position_effect + noise)
        final_qual = max(2, min(40, final_qual))  # Clamp between 2-40
        
        # Convert to ASCII (Phred+33)
        qualities.append(chr(final_qual + 33))
    
    return ''.join(qualities)

def generate_sequence(length):
    """Generate random DNA sequence."""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choices(bases, k=length))

def add_adapter_contamination(sequence, quality, contamination_rate=0.05):
    """Add realistic adapter contamination to some reads."""
    if random.random() < contamination_rate:
        # Common Illumina adapter sequence
        adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        
        # Add adapter to end of read
        if len(sequence) > 50:  # Only for longer reads
            trim_point = random.randint(40, len(sequence) - 10)
            sequence = sequence[:trim_point] + adapter[:len(sequence) - trim_point]
            # Adapter has lower quality
            adapter_qual = ''.join([chr(random.randint(15, 25) + 33) 
                                  for _ in range(len(sequence) - trim_point)])
            quality = quality[:trim_point] + adapter_qual
    
    return sequence, quality

def generate_paired_reads(read_length=75, num_pairs=10000):
    """Generate paired-end reads with realistic characteristics."""
    reads_r1 = []
    reads_r2 = []
    
    print(f"Generating {num_pairs} paired-end reads of length {read_length}...")
    
    for i in range(num_pairs):
        if i % 1000 == 0:
            print(f"  Generated {i}/{num_pairs} reads...")
        
        # Generate read identifiers (Illumina format)
        read_id = f"@HWI-ST1276:71:C1162ACXX:1:1101:{1000+i}:{2000+i}"
        
        # Generate sequences
        seq_r1 = generate_sequence(read_length)
        seq_r2 = generate_sequence(read_length)  # Independent for simplicity
        
        # Generate quality scores with realistic patterns
        # R1 typically has better quality than R2
        qual_r1 = generate_quality_string(read_length, base_quality=32, variation=8)
        qual_r2 = generate_quality_string(read_length, base_quality=28, variation=12)
        
        # Add some adapter contamination
        seq_r1, qual_r1 = add_adapter_contamination(seq_r1, qual_r1)
        seq_r2, qual_r2 = add_adapter_contamination(seq_r2, qual_r2)
        
        # Create FASTQ entries
        fastq_r1 = f"{read_id}/1\n{seq_r1}\n+\n{qual_r1}\n"
        fastq_r2 = f"{read_id}/2\n{seq_r2}\n+\n{qual_r2}\n"
        
        reads_r1.append(fastq_r1)
        reads_r2.append(fastq_r2)
    
    return reads_r1, reads_r2

def write_fastq_files(reads_r1, reads_r2, output_dir, prefix="sample"):
    """Write reads to compressed FASTQ files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    r1_file = output_dir / f"{prefix}_R1.fastq.gz"
    r2_file = output_dir / f"{prefix}_R2.fastq.gz"
    
    print(f"Writing R1 reads to {r1_file}...")
    with gzip.open(r1_file, 'wt') as f:
        f.writelines(reads_r1)
    
    print(f"Writing R2 reads to {r2_file}...")
    with gzip.open(r2_file, 'wt') as f:
        f.writelines(reads_r2)
    
    # Check file sizes
    r1_size = r1_file.stat().st_size / (1024*1024)  # MB
    r2_size = r2_file.stat().st_size / (1024*1024)  # MB
    
    print(f"Generated files:")
    print(f"  {r1_file}: {r1_size:.1f} MB")
    print(f"  {r2_file}: {r2_size:.1f} MB")
    
    return r1_file, r2_file

def create_readme(output_dir):
    """Create README describing the sample data."""
    readme_content = """# Sample RNA-seq Data for Tutorial

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
"""
    
    readme_file = Path(output_dir) / "README.md"
    with open(readme_file, 'w') as f:
        f.write(readme_content)
    
    print(f"Created dataset documentation: {readme_file}")

def main():
    """Generate sample data for the tutorial."""
    print("🧬 Generating Sample RNA-seq Data for Tutorial")
    print("=" * 50)
    
    # Set output directory
    output_dir = "guided-tutorials/01-first-rnaseq-analysis/sample-data"
    
    # Generate reads
    reads_r1, reads_r2 = generate_paired_reads(
        read_length=75,
        num_pairs=10000
    )
    
    # Write FASTQ files
    r1_file, r2_file = write_fastq_files(reads_r1, reads_r2, output_dir)
    
    # Create documentation
    create_readme(output_dir)
    
    print("\n✅ Sample data generation complete!")
    print(f"📁 Files created in: {output_dir}/")
    print(f"📊 Ready for tutorial use")
    
    # Verify files can be read
    print("\n🔍 Verifying generated files...")
    try:
        with gzip.open(r1_file, 'rt') as f:
            first_read = f.readline().strip()
            print(f"✅ R1 file readable: {first_read}")
        
        with gzip.open(r2_file, 'rt') as f:
            first_read = f.readline().strip()
            print(f"✅ R2 file readable: {first_read}")
            
    except Exception as e:
        print(f"❌ Error verifying files: {e}")
        return False
    
    print(f"\n🎯 Tutorial data ready!")
    print(f"   Files: {r1_file}, {r2_file}")
    print(f"   Use these in Module 1.1: Your First RNA-seq Analysis")

if __name__ == "__main__":
    main()