Perform STAR alignment for RNA-seq data following best practices:

1. **Check for STAR genome index** in data/reference/ directory
   - If missing, generate index using reference genome and GTF annotation
   - Use appropriate parameters for read length and species

2. **Run STAR alignment** on all paired-end FASTQ files:
   - Use 2-pass mapping for novel junction discovery
   - Output sorted BAM files
   - Generate gene counts during alignment
   - Include appropriate parameters for RNA-seq data

3. **Post-alignment processing**:
   - Index all BAM files using samtools
   - Generate alignment statistics with samtools flagstat
   - Create alignment summary report

4. **Quality control checks**:
   - Check alignment rates (should be >85% for well-annotated genomes)
   - Verify unique mapping rates (should be >70%)
   - Assess insert size distribution for paired-end data
   - Check for strand specificity if using stranded protocols

5. **Output organization**:
   - Save BAM files to results/alignments/
   - Save gene counts to results/counts/
   - Save alignment statistics to results/qc/alignment/
   - Create alignment summary report

6. **Generate alignment QC plots**:
   - Alignment rate by sample
   - Insert size distribution
   - Read mapping categories (unique, multi-mapped, unmapped)

Use optimized STAR parameters for RNA-seq:
- Enable gene counting during alignment
- Use splice junction database for better alignment
- Set appropriate multimap and mismatch parameters
- Handle paired-end data correctly

Provide clear summary of alignment results and flag any samples with concerning metrics.