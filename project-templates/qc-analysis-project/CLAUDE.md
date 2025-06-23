# Quality Control Analysis Project

## Project Overview
Comprehensive quality control workflow for bioinformatics data, including sequencing data assessment, contamination screening, and analysis validation protocols.

## Data Information
- **Data Type**: [RNA-seq/WGS/WES/ChIP-seq/ATAC-seq/other]
- **Species**: [e.g., Homo sapiens]
- **Sequencing Platform**: [Illumina/PacBio/Nanopore]
- **Read Type**: [Single-end/Paired-end]
- **Read Length**: [e.g., 150bp]
- **Number of Samples**: [Sample count]
- **Expected Coverage**: [Depth or read count]

## Quality Control Workflow

### 1. Pre-Analysis QC
- **Raw data assessment**: File integrity and format validation
- **Metadata validation**: Sample information completeness
- **Storage assessment**: Disk space and backup verification
- **Tool availability**: Required software installation check

### 2. Sequencing Data QC
- **Per-base quality scores**: FastQC analysis
- **Sequence composition**: GC content and bias assessment
- **Adapter contamination**: Trimming recommendations
- **Duplicate assessment**: PCR and optical duplicates
- **Length distribution**: Read length consistency
- **Complexity analysis**: Sequence diversity metrics

### 3. Contamination Screening
- **Species verification**: Taxonomic classification
- **Cross-contamination**: Sample mix-up detection
- **Vector contamination**: Cloning vector sequences
- **Adapter sequences**: Sequencing adapter presence
- **PhiX contamination**: Control spike-in assessment

### 4. Alignment QC (if applicable)
- **Mapping statistics**: Alignment rates and quality
- **Coverage analysis**: Depth and uniformity
- **Insert size distribution**: Paired-end metrics
- **Strand bias**: Directional sequencing protocols
- **Reference genome validation**: Version and completeness

### 5. Analysis-Specific QC
- **Expression QC**: Gene detection and correlation (RNA-seq)
- **Variant QC**: Ti/Tv ratios, allele balance (DNA-seq)
- **Peak QC**: Signal-to-noise, enrichment (ChIP-seq)
- **Fragment distribution**: Size profiles (ATAC-seq)

### 6. Report Generation
- **Summary statistics**: Key metrics across samples
- **Pass/fail assessment**: Quality threshold evaluation
- **Recommendations**: Data processing suggestions
- **Flagged samples**: Quality issues identification

## Essential Commands for Claude Code

**Context is automatically loaded!** When you start Claude Code in this project directory, the bioinformatics context documents are automatically available.

When starting QC sessions, simply describe your requirements:
```
I'm performing quality control analysis on [data type] sequencing data.

QC Requirements:
- Data type: [your data type]
- Quality standards: [your thresholds] 
- Expected metrics: [coverage, quality, etc.]
- Critical factors: [study-specific considerations]
```

The bioinformatics context (file formats, quality standards, QC tools, best practices) is already loaded from your global configuration.

## Custom Commands Available

### `/qc:fastq`
Comprehensive FASTQ quality assessment with FastQC and MultiQC

### `/qc:contamination`
Screen for various types of contamination (species, vector, adapter)

### `/qc:alignment`
Assess alignment quality and coverage metrics

### `/qc:summary`
Generate comprehensive QC report across all samples

### `/qc:outliers`
Identify samples with quality issues requiring attention

### `/qc:recommendations`
Provide data processing recommendations based on QC results

## Directory Structure
```
qc-analysis-project/
├── data/
│   ├── raw/                    # Original data files
│   ├── test/                   # Small test datasets
│   └── metadata/               # Sample information files
├── scripts/                    # QC scripts and workflows
├── results/
│   ├── fastqc/                # Per-sample FastQC outputs
│   ├── multiqc/               # Aggregated QC reports
│   ├── contamination/         # Contamination screening results
│   ├── alignment_qc/          # Alignment quality metrics
│   ├── plots/                 # QC visualization plots
│   └── reports/               # Summary reports and recommendations
├── reference/                  # Reference files for QC
├── docs/                      # QC protocols and documentation
└── .claude/
    ├── settings.json          # Project-specific settings
    └── commands/              # Custom QC commands
```

## Quality Thresholds by Data Type

### RNA-seq Quality Standards
- **Per-base quality**: Q30+ for 80%+ of bases
- **Sequence duplication**: <30% (higher acceptable for RNA-seq)
- **Adapter contamination**: <5%
- **rRNA contamination**: <10%
- **Alignment rate**: >85% for well-annotated genomes
- **Gene body coverage**: Uniform 5' to 3' distribution
- **Strand specificity**: >80% for stranded protocols

### DNA-seq (WGS/WES) Quality Standards
- **Per-base quality**: Q30+ for 80%+ of bases
- **Sequence duplication**: <20%
- **Adapter contamination**: <5%
- **Alignment rate**: >95% for WGS, >98% for WES
- **Coverage uniformity**: CV < 0.2 for WES
- **Insert size**: Consistent with library preparation

### ChIP-seq Quality Standards
- **Library complexity**: >80% of reads in peaks
- **Signal-to-noise**: >5 fold enrichment
- **Peak calling**: >1000 reproducible peaks
- **Fragment length**: Consistent with expected size
- **Cross-correlation**: Clear enrichment signal

### ATAC-seq Quality Standards
- **Fragment distribution**: Clear nucleosome patterns
- **TSS enrichment**: >5 fold at promoters
- **Library complexity**: <20% duplication
- **Peak overlap**: >70% reproducibility between replicates
- **Mitochondrial reads**: <20% of total

## Contamination Detection Protocols

### Species Verification
```bash
# Kraken2 taxonomic classification
kraken2 --db standard --threads 8 --report species_report.txt \
    --paired reads_1.fq reads_2.fq > species_classification.txt

# Expected: >95% reads classified as target species
```

### Cross-Contamination Detection
```bash
# Sample fingerprinting with common SNPs
# Compare genotype concordance between samples
# Flag samples with >5% shared variants (unexpected)
```

### Vector/Adapter Screening
```bash
# FastQ Screen against contaminant databases
fastq_screen --aligner bwa --conf fastq_screen.conf \
    --threads 8 reads.fq

# Expected: <1% reads mapping to vectors/adapters
```

## QC Visualization Standards

### Essential Plots
- **Quality score distributions**: Per-base and per-sequence
- **GC content distribution**: Compare to expected
- **Duplication levels**: Assess library complexity
- **Adapter content**: Contamination levels
- **Length distribution**: Read length consistency
- **Coverage plots**: Depth and uniformity (if aligned)

### Comparative Analysis
- **Sample correlation**: Expression or coverage correlation
- **PCA plots**: Sample clustering and outlier detection
- **Batch effect assessment**: Technical vs biological variation
- **Replicate concordance**: Technical and biological replicates

## Sample Assessment Criteria

### Pass/Fail Thresholds
**Pass**: Meets all quality thresholds for intended analysis
**Warning**: Minor quality issues, proceed with caution
**Fail**: Significant quality problems, requires reprocessing or exclusion

### Quality Categories
1. **Excellent**: Above recommended thresholds
2. **Good**: Meets minimum requirements
3. **Marginal**: Below optimal but potentially usable
4. **Poor**: Requires reprocessing or exclusion

## Common QC Issues and Solutions

### Low Quality Scores
- **Cause**: Sequencing chemistry or cycle issues
- **Solution**: Quality trimming, consider re-sequencing
- **Impact**: Reduced accuracy, increased false positives

### High Duplication Rates
- **Cause**: Over-amplification, low input material
- **Solution**: Assess library complexity, mark duplicates
- **Impact**: Reduced effective coverage, biased representation

### Adapter Contamination
- **Cause**: Insert size shorter than read length
- **Solution**: Adapter trimming with cutadapt/trimmomatic
- **Impact**: Alignment artifacts, false variants

### Uneven Coverage
- **Cause**: GC bias, capture efficiency issues
- **Solution**: Bias correction, coverage normalization
- **Impact**: Missed regions, variable sensitivity

### Batch Effects
- **Cause**: Technical differences between runs
- **Solution**: Batch correction, balanced design
- **Impact**: Confounded biological signals

## Automated QC Pipeline

### Daily QC Checks
```bash
# Automated daily QC for new data
for sample in new_samples/*; do
    fastqc $sample
    fastq_screen $sample
done
multiqc results/
```

### Weekly QC Review
- Aggregate QC metrics across all samples
- Identify trending quality issues
- Update QC thresholds based on data
- Generate lab QC summary report

### Monthly QC Audit
- Review QC protocols and thresholds
- Assess tool performance and updates
- Update reference databases
- Train team on new QC methods

## Reporting Templates

### Sample-Level Report
```
Sample: [sample_id]
Status: [PASS/WARNING/FAIL]
Quality Score: [mean Q score]
Duplication Rate: [percentage]
Contamination: [percentage]
Recommendations: [specific actions needed]
```

### Project-Level Summary
```
Project: [project_name]
Total Samples: [count]
Pass Rate: [percentage]
Common Issues: [list top 3]
Overall Recommendation: [proceed/reprocess/investigate]
```

## Important Notes
- Run QC before any downstream analysis
- Maintain QC logs for all projects
- Update QC thresholds based on study requirements
- Document any QC protocol deviations
- Archive QC reports for reproducibility

## Team Contacts
- **QC Manager**: [Name and contact]
- **Bioinformatics Lead**: [Name and contact]
- **Sequencing Core**: [Core facility contact]
- **Technical Support**: [IT/computational support]

---

*This CLAUDE.md file provides comprehensive QC protocols that should be adapted based on specific data types and study requirements.*