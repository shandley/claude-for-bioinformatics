# Variant Calling Project

## Project Overview
Comprehensive variant calling pipeline following GATK best practices for whole genome sequencing (WGS) or whole exome sequencing (WES) data.

## Data Information
- **Species**: [e.g., Homo sapiens]
- **Reference Genome**: [e.g., GRCh38]
- **Data Type**: [WGS/WES] paired-end sequencing
- **Read Length**: [e.g., 150bp]
- **Samples**: [Number of samples and study design]
- **Population**: [Ancestry/population information]
- **Coverage**: [Expected coverage depth]

## Analysis Workflow

### 1. Quality Control
- FastQC analysis of raw reads
- MultiQC aggregated report
- Contamination screening
- Coverage assessment

### 2. Read Preprocessing
- Adapter trimming (if needed)
- Quality filtering
- Duplicate marking preparation

### 3. Alignment
- BWA-MEM alignment to reference genome
- Coordinate sorting and indexing
- Alignment quality assessment

### 4. Post-Alignment Processing
- Mark duplicate reads with Picard
- Base quality score recalibration (BQSR)
- Alignment refinement around indels

### 5. Variant Calling
- HaplotypeCaller in GVCF mode
- Joint genotyping with GenotypeGVCFs
- Variant quality assessment

### 6. Variant Filtering
- Hard filtering or VQSR based on cohort size
- Population-specific filters
- Quality metric thresholds

### 7. Annotation and Interpretation
- VEP or SnpEff annotation
- Population frequency annotation
- Clinical significance annotation
- Functional impact prediction

## Essential Commands for Claude Code

**Context is automatically loaded!** When you start Claude Code in this project directory, the bioinformatics context documents are automatically available.

When starting analysis sessions, simply describe your goal:
```
I'm working on variant calling analysis following GATK best practices.

Project details:
- Species: [your species]
- Reference: [your reference genome version] 
- Data type: [WGS/WES]
- Samples: [your sample information]
- Goal: [germline/somatic variant detection]
```

The bioinformatics context (file formats, tools, GATK best practices, quality standards) is already loaded from your global configuration.

## Custom Commands Available

### `/variants:qc`
Run comprehensive quality control on FASTQ files

### `/variants:align`
Perform BWA-MEM alignment with GATK-compatible settings

### `/variants:preprocess`
Mark duplicates and perform base quality recalibration

### `/variants:call`
Run GATK HaplotypeCaller in GVCF mode

### `/variants:joint-genotype`
Perform joint genotyping across all samples

### `/variants:filter`
Apply hard filters or VQSR based on cohort size

### `/variants:annotate`
Annotate variants with functional and population information

## Directory Structure
```
variant-calling-project/
├── data/
│   ├── raw/                    # Original FASTQ files
│   ├── processed/              # Quality-controlled reads
│   └── reference/              # Reference genome, indices, known sites
├── scripts/                    # Analysis scripts and workflows
├── results/
│   ├── qc/                    # Quality control reports
│   ├── alignments/            # BAM files and indices
│   ├── gvcfs/                 # Individual sample GVCFs
│   ├── variants/              # Joint-called and filtered VCFs
│   ├── annotations/           # Annotated variant files
│   └── reports/               # Analysis summaries and plots
├── docs/                      # Analysis documentation
└── .claude/
    ├── settings.json          # Project-specific settings
    └── commands/              # Custom variant calling commands
```

## Quality Control Thresholds

### FASTQ Quality Standards
- **Per-base quality**: Q30+ for 80%+ of bases
- **Sequence duplication**: <20%
- **Adapter contamination**: <5%
- **Insert size**: Appropriate for sequencing protocol

### Alignment Quality Standards
- **Alignment rate**: >95% for WGS, >98% for WES
- **Properly paired**: >95% for paired-end data
- **Coverage uniformity**: CV < 0.2 for WES
- **Mean coverage**: 30x for WGS, 100x for WES

### Variant Quality Standards
- **QUAL score**: >30 for high-confidence variants
- **Depth (DP)**: 10-100x for WGS, 20-200x for WES
- **Genotype Quality (GQ)**: >20
- **Allele Balance**: 0.2-0.8 for heterozygous calls

## Reference Files Required
- **Reference genome FASTA**: [e.g., GRCh38.primary_assembly.genome.fa]
- **BWA index**: Generated from reference genome
- **Sequence dictionary**: Created with Picard
- **Known sites VCFs**:
  - dbSNP (latest version)
  - Mills and 1000G gold standard indels
  - 1000 Genomes Phase 3 SNPs and indels

## Sample Information Template
```
sample_id,population,sex,phenotype,batch,notes
sample_01,EUR,F,control,batch1,
sample_02,EUR,M,control,batch1,
sample_03,AFR,F,case,batch2,
sample_04,ASN,M,case,batch2,
```

## GATK Best Practices Parameters

### BWA-MEM Alignment
```bash
bwa mem -t [threads] -M -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
    reference.fa read1.fq read2.fq | \
    samtools sort -@ [threads] -o sample.bam
```

### Mark Duplicates
```bash
gatk MarkDuplicates \
    -I input.bam \
    -O marked_duplicates.bam \
    -M metrics.txt
```

### Base Quality Score Recalibration
```bash
# Generate recalibration table
gatk BaseRecalibrator \
    -I marked_duplicates.bam \
    -R reference.fa \
    --known-sites dbsnp.vcf \
    --known-sites mills.vcf \
    -O recal_data.table

# Apply recalibration
gatk ApplyBQSR \
    -I marked_duplicates.bam \
    -R reference.fa \
    --bqsr-recal-file recal_data.table \
    -O recalibrated.bam
```

### Variant Calling
```bash
gatk HaplotypeCaller \
    -I recalibrated.bam \
    -R reference.fa \
    -ERC GVCF \
    -O sample.g.vcf.gz
```

### Joint Genotyping
```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V cohort.g.vcf.gz \
    -O cohort.vcf.gz
```

## Filtering Strategies

### Hard Filters (Small Cohorts)
**SNPs**:
- QD < 2.0
- QUAL < 30.0
- SOR > 3.0
- FS > 60.0
- MQ < 40.0
- MQRankSum < -12.5
- ReadPosRankSum < -8.0

**Indels**:
- QD < 2.0
- QUAL < 30.0
- FS > 200.0
- ReadPosRankSum < -20.0

### VQSR (Large Cohorts)
- Use when >30 exomes or >10 genomes
- Truth/training sets: HapMap, Omni, 1000G, dbSNP
- Sensitivity thresholds: 99.7% for SNPs, 99.0% for indels

## Common Troubleshooting

### Low Alignment Rates
- Check FASTQ quality and adapter contamination
- Verify reference genome matches sample ancestry
- Assess for sample contamination

### High Duplicate Rates
- Check for over-amplification during library prep
- Consider optical duplicate marking
- Assess impact on coverage uniformity

### Poor Variant Quality
- Check base quality recalibration effectiveness
- Assess coverage depth and uniformity
- Verify known sites files are appropriate

### Batch Effects
- Check for systematic differences between batches
- Consider batch as covariate in analysis
- Use appropriate population controls

## Population Considerations

### Ancestry Matching
- Use appropriate reference panels
- Consider population-specific allele frequencies
- Account for admixed populations

### Known Sites Selection
- Use population-matched training sets when available
- Include relevant population databases
- Consider local variant databases

## Analysis Timeline
- **Day 1-2**: Quality control and alignment
- **Day 3-4**: Preprocessing and BQSR
- **Day 5-6**: Variant calling and joint genotyping
- **Day 7**: Filtering and quality assessment
- **Day 8-9**: Annotation and interpretation
- **Day 10**: Report generation and validation

## Important Notes
- Always validate key findings with orthogonal methods
- Maintain detailed logs of all parameters used
- Check for sex chromosome consistency
- Monitor computational resource usage
- Back up critical intermediate files

## Team Contacts
- **Principal Investigator**: [Name and contact]
- **Bioinformatics Lead**: [Name and contact]
- **Clinical Geneticist**: [Name and contact for clinical interpretation]
- **Technical Support**: [Core facility or IT contact]

---

*This CLAUDE.md file should be customized for each specific variant calling project while maintaining adherence to GATK best practices.*