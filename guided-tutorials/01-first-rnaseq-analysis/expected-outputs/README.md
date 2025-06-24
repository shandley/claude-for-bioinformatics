# Expected Outputs for Module 1.1

## What You Should See After Completing the Tutorial

After successfully running the FastQC and MultiQC analysis on the sample data, you should have generated the following files and results:

### File Structure
```
results/qc/
├── sample_R1_fastqc.html      # FastQC report for R1 reads
├── sample_R1_fastqc.zip       # FastQC data files for R1
├── sample_R2_fastqc.html      # FastQC report for R2 reads  
├── sample_R2_fastqc.zip       # FastQC data files for R2
├── multiqc_report.html        # Combined MultiQC report
└── multiqc_data/              # MultiQC data directory
    ├── multiqc.log            # MultiQC processing log
    ├── multiqc_data.json      # Machine-readable data
    ├── multiqc_fastqc.txt     # FastQC summary data
    └── multiqc_general_stats.txt  # General statistics
```

## Expected Quality Control Results

### Overall Assessment
✅ **PASS**: This sample data should produce mostly "PASS" results with some educational "WARN" flags

### FastQC Module Results

#### ✅ Basic Statistics (PASS)
- **Total Sequences**: 10,000 per file
- **Sequence Length**: 75 bp
- **%GC**: ~50% (normal for synthetic data)
- **Total Duplicated Reads**: ~15-25% (normal for RNA-seq)

#### ✅ Per-base Sequence Quality (PASS)
- Quality scores should be mostly >28 (green zone)
- Slight quality decline toward 3' end (typical Illumina pattern)
- R1 should have better quality than R2

#### ⚠️ Per-sequence Quality Scores (WARN)
- Some reads with lower quality scores (educational value)
- Demonstrates quality distribution interpretation

#### ✅ Per-base Sequence Content (PASS/WARN)
- May show slight bias in first few bases (normal)
- Good example for discussing library prep artifacts

#### ✅ Per-sequence GC Content (PASS)
- Normal distribution around 50%
- Good for discussing species-specific expectations

#### ⚠️ Per-base N Content (PASS)
- Should be very low (<1%)
- Good example of high-quality data

#### ⚠️ Sequence Length Distribution (WARN)
- All reads 75bp (uniform length)
- May trigger warning for lack of variation

#### ⚠️ Sequence Duplication Levels (WARN)
- Moderate duplication (15-25%)
- Educational discussion about RNA-seq duplication patterns

#### ⚠️ Overrepresented Sequences (WARN)
- May show some adapter sequences
- Good for discussing adapter trimming needs

#### ⚠️ Adapter Content (WARN)
- ~5% adapter contamination at read ends
- Perfect for learning about adapter detection and trimming

### MultiQC Summary Statistics

Expected values in the MultiQC general statistics table:

| Sample    | Total Sequences | GC% | Avg Length | % Duplicates | Adapter% |
|-----------|----------------|-----|------------|-------------|----------|
| sample_R1 | 10,000         | ~50 | 75         | 15-25       | 3-7      |
| sample_R2 | 10,000         | ~50 | 75         | 15-25       | 3-7      |

### Key Learning Points

#### 1. Quality Score Interpretation
- **Q30**: 99.9% accuracy (1 in 1000 error rate)
- **Q20**: 99% accuracy (1 in 100 error rate)
- Most reads should be above Q20, preferably Q30+

#### 2. Typical RNA-seq Patterns
- **Duplication**: Higher than DNA-seq due to expression levels
- **GC bias**: May reflect library prep or species characteristics
- **Quality decline**: Normal toward 3' end of reads

#### 3. When to Worry
- **Major quality drops**: <Q20 for large portions of reads
- **High adapter content**: >10% suggests trimming needed
- **Extreme GC bias**: May indicate contamination
- **Very high duplication**: >50% suggests PCR over-amplification

#### 4. Next Steps Decision Making
Based on these results, students should conclude:
- ✅ **Data quality**: Generally good, suitable for analysis
- ⚠️ **Adapter trimming**: Recommended due to 5% contamination
- ✅ **Proceed to alignment**: Quality sufficient for mapping
- 📝 **Document parameters**: Record quality thresholds used

## Troubleshooting Expected Results

### If You Don't See These Results

#### Different Numbers
- **Slight variations are normal**: Random generation creates minor differences
- **Major differences**: Check if FastQC ran on correct files
- **Missing data**: Verify FastQC completed successfully

#### No MultiQC Report
- **Check FastQC outputs exist**: Should see .zip files in results/qc/
- **Run MultiQC again**: `multiqc results/qc/ -o results/qc/ -f`
- **Check MultiQC log**: Look at `multiqc_data/multiqc.log` for errors

#### Unexpected Quality Issues
- **Much worse quality**: Verify you're using the tutorial sample data
- **Much better quality**: Also check you're using correct files
- **Different patterns**: Compare with this expected output guide

## Educational Value

This sample dataset is specifically designed to:

1. **Show realistic patterns** without requiring large computational resources
2. **Include common issues** (adapter contamination, quality variation) for learning
3. **Produce interpretable results** in under 2 minutes processing time
4. **Demonstrate decision-making** about data quality and next steps
5. **Practice troubleshooting** with realistic but manageable complexity

The results should give students confidence in:
- Reading and interpreting QC reports
- Understanding what "good enough" data looks like
- Making informed decisions about data processing steps
- Recognizing when additional quality control steps are needed

---

*If your results differ significantly from these expectations, please check the troubleshooting guide or ask for help in the GitHub discussions.*