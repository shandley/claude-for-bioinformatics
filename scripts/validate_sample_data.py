#!/usr/bin/env python3
"""
Validate the generated sample data for tutorial use.

This script performs basic validation checks on the sample FASTQ files
to ensure they will work correctly in the tutorial environment.
"""

import gzip
import os
from pathlib import Path

def validate_fastq_format(filepath):
    """Validate basic FASTQ format requirements."""
    print(f"📁 Validating {filepath}")
    
    try:
        with gzip.open(filepath, 'rt') as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                if i >= 7:  # Read first 2 records (8 lines)
                    break
            
            # Check we have enough lines
            if len(lines) < 8:
                print(f"❌ File too short: only {len(lines)} lines")
                return False
            
            # Validate first record
            if not lines[0].startswith('@'):
                print(f"❌ Invalid header line: {lines[0][:50]}")
                return False
            
            if not lines[2].startswith('+'):
                print(f"❌ Invalid separator line: {lines[2][:50]}")
                return False
            
            if len(lines[1]) != len(lines[3]):
                print(f"❌ Sequence/quality length mismatch: {len(lines[1])} vs {len(lines[3])}")
                return False
            
            # Check quality scores are in valid range
            for qual_char in lines[3]:
                qual_score = ord(qual_char) - 33
                if qual_score < 0 or qual_score > 93:
                    print(f"❌ Invalid quality score: {qual_char} ({qual_score})")
                    return False
            
            print(f"✅ Format validation passed")
            print(f"   Read length: {len(lines[1])} bp")
            print(f"   Quality range: {min(ord(c)-33 for c in lines[3])}-{max(ord(c)-33 for c in lines[3])}")
            
            return True
            
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return False

def count_reads(filepath):
    """Count total number of reads in FASTQ file."""
    try:
        with gzip.open(filepath, 'rt') as f:
            count = 0
            line_num = 0
            for line in f:
                if line_num % 4 == 0 and line.startswith('@'):
                    count += 1
                line_num += 1
        print(f"📊 Read count: {count:,}")
        return count
    except Exception as e:
        print(f"❌ Error counting reads: {e}")
        return 0

def check_file_size(filepath):
    """Check file size is appropriate for tutorial."""
    size_mb = filepath.stat().st_size / (1024 * 1024)
    print(f"📏 File size: {size_mb:.1f} MB")
    
    if size_mb > 50:
        print(f"⚠️  File might be too large for quick tutorial processing")
        return False
    elif size_mb < 0.1:
        print(f"⚠️  File might be too small for meaningful analysis")
        return False
    else:
        print(f"✅ File size appropriate for tutorial")
        return True

def validate_paired_reads(r1_file, r2_file):
    """Validate that R1 and R2 files are properly paired."""
    print(f"🔗 Validating read pairing...")
    
    try:
        with gzip.open(r1_file, 'rt') as f1, gzip.open(r2_file, 'rt') as f2:
            for i in range(5):  # Check first 5 read pairs
                r1_header = f1.readline().strip()
                r1_seq = f1.readline().strip()
                r1_plus = f1.readline().strip()
                r1_qual = f1.readline().strip()
                
                r2_header = f2.readline().strip()
                r2_seq = f2.readline().strip()
                r2_plus = f2.readline().strip()
                r2_qual = f2.readline().strip()
                
                # Extract read identifiers (remove /1 and /2)
                r1_id = r1_header.split('/')[0] if '/' in r1_header else r1_header.split()[0]
                r2_id = r2_header.split('/')[0] if '/' in r2_header else r2_header.split()[0]
                
                if r1_id != r2_id:
                    print(f"❌ Read pairing mismatch at pair {i+1}")
                    print(f"   R1: {r1_header}")
                    print(f"   R2: {r2_header}")
                    return False
        
        print(f"✅ Read pairing validation passed")
        return True
        
    except Exception as e:
        print(f"❌ Error validating pairing: {e}")
        return False

def main():
    """Run all validation checks."""
    print("🧬 Validating Sample RNA-seq Data")
    print("=" * 50)
    
    # File paths
    sample_dir = Path("guided-tutorials/01-first-rnaseq-analysis/sample-data")
    r1_file = sample_dir / "sample_R1.fastq.gz"
    r2_file = sample_dir / "sample_R2.fastq.gz"
    
    # Check files exist
    if not r1_file.exists():
        print(f"❌ R1 file not found: {r1_file}")
        return False
    
    if not r2_file.exists():
        print(f"❌ R2 file not found: {r2_file}")
        return False
    
    print(f"📂 Found sample data files:")
    print(f"   {r1_file}")
    print(f"   {r2_file}")
    print()
    
    # Validation checks
    all_passed = True
    
    # Validate R1 file
    print("🔍 Validating R1 file...")
    if not validate_fastq_format(r1_file):
        all_passed = False
    if not check_file_size(r1_file):
        all_passed = False
    r1_count = count_reads(r1_file)
    print()
    
    # Validate R2 file
    print("🔍 Validating R2 file...")
    if not validate_fastq_format(r2_file):
        all_passed = False
    if not check_file_size(r2_file):
        all_passed = False
    r2_count = count_reads(r2_file)
    print()
    
    # Check read counts match
    if r1_count != r2_count:
        print(f"❌ Read count mismatch: R1={r1_count}, R2={r2_count}")
        all_passed = False
    else:
        print(f"✅ Read counts match: {r1_count:,} read pairs")
    
    # Validate pairing
    if not validate_paired_reads(r1_file, r2_file):
        all_passed = False
    
    print()
    
    # Final assessment
    if all_passed:
        print("🎉 All validation checks passed!")
        print("📚 Sample data is ready for tutorial use")
        print()
        print("Tutorial performance estimates:")
        print(f"   FastQC processing time: ~30-60 seconds")
        print(f"   MultiQC processing time: ~5-10 seconds")
        print(f"   Total tutorial time: ~45-60 minutes")
        return True
    else:
        print("❌ Some validation checks failed")
        print("🔧 Please regenerate sample data or fix issues")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)