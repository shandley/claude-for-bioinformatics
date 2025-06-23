# RNA-seq Analysis Project

## Project Overview
Comprehensive RNA-seq analysis pipeline for differential expression analysis, including quality control, alignment, quantification, and statistical analysis.

## Data Information
- **Species**: [e.g., Homo sapiens]
- **Reference Genome**: [e.g., GRCh38]
- **Data Type**: Paired-end RNA-seq
- **Read Length**: [e.g., 150bp]
- **Samples**: [Number of samples and experimental design]
- **Conditions**: [Treatment groups, time points, etc.]

## Analysis Workflow

### 1. Quality Control
- FastQC analysis of raw reads
- MultiQC aggregated report
- Adapter contamination assessment
- Read quality and duplication analysis

### 2. Preprocessing
- Adapter trimming (if needed)
- Quality filtering
- rRNA contamination removal (if needed)

### 3. Alignment and Quantification
- STAR alignment to reference genome
- Gene-level quantification with featureCounts
- OR Pseudo-alignment with Salmon for transcript quantification

### 4. Differential Expression Analysis
- Count normalization and filtering
- Statistical analysis with DESeq2
- Multiple testing correction
- Fold change and significance thresholds

### 5. Functional Analysis
- Gene Ontology (GO) enrichment
- KEGG pathway analysis
- Gene Set Enrichment Analysis (GSEA)

### 6. Visualization and Reporting
- PCA and sample clustering
- MA plots and volcano plots
- Heatmaps of significant genes
- Summary report generation

## Essential Commands for Claude Code

**Context is automatically loaded!** When you start Claude Code in this project directory, the bioinformatics context documents are automatically available.

When starting analysis sessions, simply describe your goal:
```
I'm working on RNA-seq differential expression analysis.

Project details:
- Species: [your species] 
- Reference: [your reference genome]
- Samples: [your sample information]
- Goal: [your research question]
```

The bioinformatics context (file formats, tools, quality standards, best practices) is already loaded from your global configuration.

## Custom Commands Available

Use these slash commands for common RNA-seq tasks:

### `/rnaseq:qc`
Run comprehensive quality control on all FASTQ files in data/raw/

### `/rnaseq:align`
Perform STAR alignment with optimized parameters for this project

### `/rnaseq:quantify`
Generate gene counts using featureCounts with project GTF

### `/rnaseq:deseq2`
Run differential expression analysis with DESeq2

### `/rnaseq:plots`
Generate standard visualization plots (PCA, volcano, MA)

## Directory Structure
```
rnaseq-project/
├── data/
│   ├── raw/                    # Original FASTQ files
│   ├── processed/              # Quality-controlled reads
│   └── reference/              # Reference genome and annotations
├── scripts/                    # Analysis scripts and workflows
├── results/
│   ├── qc/                    # Quality control reports
│   ├── alignments/            # BAM files and indices
│   ├── counts/                # Gene expression matrices
│   ├── differential/          # DE analysis results
│   ├── functional/            # GO/KEGG enrichment results
│   └── plots/                 # Visualization outputs
├── docs/                      # Analysis documentation
└── .claude/
    ├── settings.json          # Project-specific settings
    └── commands/              # Custom RNA-seq commands
```

## Quality Control Thresholds

### FASTQ Quality Standards
- **Per-base quality**: Q30+ for 80%+ of bases
- **Sequence duplication**: <20% for RNA-seq
- **GC content**: Within expected range for species
- **Adapter contamination**: <5%

### Alignment Quality Standards
- **Alignment rate**: >85% for well-annotated genomes
- **Unique mapping**: >70% of aligned reads
- **Strand specificity**: >80% for stranded protocols
- **Gene body coverage**: Uniform 5' to 3' distribution

### Expression Analysis Standards
- **Gene filtering**: CPM > 1 in at least 3 samples
- **Normalization**: DESeq2 size factor normalization
- **Significance**: padj < 0.05 (FDR corrected)
- **Effect size**: |log2FC| > 1.0 for biological significance

## Reference Files Required
- **Reference genome FASTA**: [e.g., GRCh38.primary_assembly.genome.fa]
- **Gene annotation GTF**: [e.g., gencode.v44.primary_assembly.annotation.gtf]
- **STAR genome index**: Generated from reference files
- **Gene ID mapping**: For functional analysis

## Sample Metadata Template
```
sample_id,condition,replicate,batch,notes
sample_01,control,1,batch1,
sample_02,control,2,batch1,
sample_03,control,3,batch2,
sample_04,treatment,1,batch1,
sample_05,treatment,2,batch1,
sample_06,treatment,3,batch2,
```

## Common Troubleshooting

### Low Alignment Rates
- Check reference genome version matches data
- Verify read orientation (forward/reverse)
- Assess adapter contamination levels

### High Duplication Rates
- Expected for RNA-seq due to transcript abundance
- Consider biological vs technical replication
- Check for amplification artifacts

### Poor Sample Clustering
- Check batch effects in metadata
- Verify sample labeling accuracy
- Consider additional covariates

### Low Gene Detection
- Assess sequencing depth adequacy
- Check gene annotation completeness
- Verify library preparation quality

## Key Parameters for This Project

### STAR Alignment
```bash
--readFilesCommand zcat
--outSAMtype BAM SortedByCoordinate
--quantMode GeneCounts
--outFilterType BySJout
--outFilterMultimapNmax 20
--alignSJoverhangMin 8
--alignSJDBoverhangMin 1
--outFilterMismatchNmax 999
--outFilterMismatchNoverReadLmax 0.04
--alignIntronMin 20
--alignIntronMax 1000000
--alignMatesGapMax 1000000
```

### featureCounts
```bash
-p -B -C -t exon -g gene_id
-T [number of threads]
-a [GTF file]
-o [output file]
[BAM files]
```

### DESeq2 Design
```r
# Adjust based on experimental design
design = ~ condition
# Or with batch effects:
design = ~ batch + condition
```

## Analysis Timeline
- **Day 1**: Quality control and preprocessing
- **Day 2**: Alignment and quantification  
- **Day 3**: Differential expression analysis
- **Day 4**: Functional analysis and visualization
- **Day 5**: Report generation and validation

## Important Notes
- Always check intermediate outputs for quality
- Maintain analysis log with parameters used
- Back up critical results before major changes
- Document any deviations from standard protocol
- Validate key findings with independent methods

## Team Contacts
- **Principal Investigator**: [Name and contact]
- **Bioinformatics Lead**: [Name and contact]
- **Technical Support**: [Core facility or IT contact]

---

*This CLAUDE.md file should be customized for each specific RNA-seq project while maintaining the core structure and quality standards.*