Run comprehensive quality control analysis on RNA-seq data:

1. **Find all FASTQ files** in data/raw/ directory
2. **Run FastQC** on each file with appropriate parameters
3. **Generate MultiQC report** aggregating all QC results
4. **Check for common issues**:
   - Adapter contamination levels
   - Per-base quality scores
   - Sequence duplication rates
   - GC content distribution
   - Read length distribution
5. **Create QC summary** with pass/fail status for each sample
6. **Flag problematic samples** that may need additional preprocessing
7. **Save all outputs** to results/qc/ directory with timestamps
8. **Generate recommendations** for next steps based on QC results

Use RNA-seq specific quality thresholds:
- Per-base quality: Q30+ for 80%+ of bases
- Sequence duplication: <20% (higher acceptable for RNA-seq)
- Adapter contamination: <5%
- Check strand specificity if using stranded protocols

Provide clear interpretation of results and suggest preprocessing steps if needed.