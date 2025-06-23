# Bioinformatics Context Reference Guide

## Table of Contents
1. [File Formats](#file-formats)
2. [Software Tools](#software-tools)
3. [Databases](#databases)
4. [Reference Genomes & Standards](#reference-genomes--standards)
5. [Programming Languages & Libraries](#programming-languages--libraries)
6. [Quality Control Standards](#quality-control-standards)
7. [Common Commands & One-liners](#common-commands--one-liners)
8. [File Extensions Quick Reference](#file-extensions-quick-reference)
9. [Hardware & Resource Requirements](#hardware--resource-requirements)
10. [Best Practices & Standards](#best-practices--standards)

---

## File Formats

### Sequence Data Formats

#### FASTA (.fasta, .fa, .fas, .fna, .ffn, .faa, .frn)
- **Purpose:** Store nucleotide or protein sequences
- **Structure:** Header line (starts with `>`) + sequence data
- **Usage:** Reference genomes, protein sequences, database submissions
- **Tools:** Compatible with most bioinformatics software

#### FASTQ (.fastq, .fq)
- **Purpose:** Store raw sequencing reads with quality scores
- **Structure:** 4 lines per read (header, sequence, +, quality)
- **Quality encoding:** Phred+33 (Sanger), Phred+64 (Illumina 1.3+)
- **Compression:** Usually gzipped (.fastq.gz, .fq.gz)

### Alignment Formats

#### SAM (.sam) / BAM (.bam)
- **SAM:** Sequence Alignment/Map (text format)
- **BAM:** Binary compressed version of SAM
- **Structure:** Header + alignment records (11 mandatory fields)
- **Index:** BAI files (.bam.bai) for rapid access
- **Tools:** samtools, Picard, GATK

#### CRAM (.cram)
- **Purpose:** Reference-based compression of alignment data
- **Advantage:** Smaller file size than BAM
- **Usage:** Long-term storage, newer format

### Variant Formats

#### VCF (.vcf) / BCF (.bcf)
- **VCF:** Variant Call Format (text)
- **BCF:** Binary compressed VCF
- **Purpose:** Store genetic variants (SNPs, indels, structural variants)
- **Index:** Tabix (.tbi) or CSI (.csi) files
- **Compression:** Usually bgzipped (.vcf.gz)

### Annotation Formats

#### GFF3 (.gff3, .gff)
- **Purpose:** Gene Feature Format for genomic annotations
- **Structure:** Tab-delimited with 9 columns
- **Usage:** Gene predictions, functional annotations

#### GTF (.gtf)
- **Purpose:** Gene Transfer Format (subset of GFF2)
- **Usage:** Transcript annotations, RNA-seq analysis
- **Tools:** StringTie, Cufflinks, featureCounts

#### BED (.bed)
- **Purpose:** Browser Extensible Data for genomic intervals
- **Structure:** Tab-delimited (chromosome, start, end, ...)
- **Variants:** BED3, BED6, BED12 formats
- **Usage:** Peak calling, genomic regions, annotations

### Other Important Formats

#### HDF5 (.h5, .hdf5)
- **Purpose:** Hierarchical data format for large datasets
- **Usage:** Single-cell data, imaging data, large matrices
- **Tools:** h5py, rhdf5, AnnData

#### AnnData (.h5ad)
- **Purpose:** Annotated data format for single-cell analysis
- **Structure:** HDF5-based with observations, variables, metadata
- **Tools:** scanpy, Seurat (via conversion)

#### Phylogenetic Formats
- **Newick (.nwk, .tre):** Tree topology with branch lengths
- **PHYLIP (.phy):** Multiple sequence alignment format
- **NEXUS (.nex):** Extended format with metadata
- **Stockholm (.sto):** Alignment format with secondary structure

---

## Software Tools

### Sequence Analysis

#### BLAST Suite
- **blastn:** Nucleotide-nucleotide BLAST
- **blastp:** Protein-protein BLAST
- **blastx:** Nucleotide query vs protein database
- **tblastn:** Protein query vs translated nucleotide database
- **Usage:** Sequence similarity, homolog identification

#### Multiple Sequence Alignment
- **MUSCLE:** Fast multiple sequence alignment
- **MAFFT:** Accurate alignment for large datasets
- **ClustalW/ClustalX:** Traditional MSA tools
- **T-Coffee:** Consistency-based alignment

### Read Processing & Alignment

#### Quality Control
- **FastQC:** Quality assessment of sequencing data
- **MultiQC:** Aggregate QC reports from multiple tools
- **fastp:** All-in-one preprocessing tool
- **Trimmomatic:** Read trimming and filtering

#### Read Alignment
- **BWA:** Burrows-Wheeler Aligner (DNA-seq)
  - `bwa mem`: Main algorithm for modern reads
  - `bwa aln`: For shorter reads
- **STAR:** Spliced Transcripts Alignment to Reference (RNA-seq)
- **HISAT2:** Hierarchical indexing for RNA-seq
- **Bowtie2:** Fast alignment for DNA-seq
- **minimap2:** Long-read and splice-aware alignment

### Variant Calling

#### GATK (Genome Analysis Toolkit)
- **HaplotypeCaller:** Variant calling
- **GenotypeGVCFs:** Joint genotyping
- **BaseRecalibrator:** Quality score recalibration
- **Best Practices:** Standardized workflows

#### Other Variant Callers
- **FreeBayes:** Haplotype-based variant detection
- **VarScan:** Somatic and germline variant calling
- **Mutect2:** Somatic variant calling in cancer

### RNA-seq Analysis

#### Quantification
- **Salmon:** Transcript quantification with quasi-mapping
- **Kallisto:** Fast transcript quantification
- **featureCounts:** Gene-level counting from alignments
- **HTSeq:** Python-based read counting

#### Differential Expression
- **DESeq2:** R package for differential expression
- **edgeR:** R package for count data analysis
- **limma:** Linear models for microarray/RNA-seq data

### Single-cell Analysis

#### Python Ecosystem
- **scanpy:** Single-cell analysis in Python
- **CellRanger:** 10X Genomics pipeline
- **Velocyto:** RNA velocity analysis
- **scVI:** Deep generative models

#### R Ecosystem
- **Seurat:** Comprehensive single-cell toolkit
- **SingleCellExperiment:** Data structure for single-cell
- **Monocle:** Trajectory analysis
- **scran:** Quality control and normalization

### Genome Assembly

#### Short-read Assemblers
- **SPAdes:** Versatile genome assembler
- **Velvet:** De Bruijn graph assembler
- **ABySS:** Large genome assembly

#### Long-read Assemblers
- **Canu:** PacBio and Nanopore assembly
- **Flye:** Fast long-read assembler
- **Hifiasm:** PacBio HiFi assembler

### Utilities

#### File Format Conversion
- **samtools:** SAM/BAM/CRAM manipulation
- **bcftools:** VCF/BCF manipulation
- **bedtools:** BED file operations
- **seqkit:** FASTA/FASTQ processing
- **bioawk:** AWK for bioinformatics formats

---

## Databases

### NCBI Resources

#### Primary Databases
- **GenBank:** Nucleotide sequence database
- **RefSeq:** Curated reference sequences
- **PubMed:** Biomedical literature database
- **dbSNP:** Single nucleotide polymorphism database
- **dbGaP:** Genotype and phenotype database
- **SRA:** Sequence Read Archive
- **Gene:** Gene-specific information portal

#### Specialized Databases
- **ClinVar:** Clinical variant database
- **GTR:** Genetic Testing Registry
- **OMIM:** Online Mendelian Inheritance in Man
- **PMC:** PubMed Central full-text articles

### European Resources (EMBL-EBI)

#### Core Databases
- **Ensembl:** Genome browser and annotation
- **UniProt:** Protein sequence and functional information
- **ChEMBL:** Bioactive compounds database
- **ArrayExpress:** Gene expression database
- **ENA:** European Nucleotide Archive

#### Specialized Resources
- **Pfam:** Protein family database
- **InterPro:** Protein domain database
- **STRING:** Protein-protein interaction networks
- **Reactome:** Pathway database

### Other Important Databases

#### Genomics
- **UCSC Genome Browser:** Interactive genome visualization
- **1000 Genomes:** Population genetics reference
- **gnomAD:** Genome aggregation database
- **COSMIC:** Catalogue of somatic mutations in cancer

#### Protein Structure
- **PDB:** Protein Data Bank (3D structures)
- **AlphaFold:** AI-predicted protein structures
- **ChimeraX:** Molecular visualization

#### Model Organisms
- **SGD:** Saccharomyces Genome Database (yeast)
- **FlyBase:** Drosophila genome database
- **WormBase:** C. elegans genome database
- **TAIR:** Arabidopsis information resource
- **MGI:** Mouse genome informatics

---

## Reference Genomes & Standards

### Human Reference Genomes

#### Current Standards
- **GRCh38/hg38:** Current human reference (2013-)
- **GRCh37/hg19:** Previous reference (still widely used)
- **T2T-CHM13:** Telomere-to-telomere complete genome (2022)

#### Key Features
- **Chromosome naming:** chr1-chr22, chrX, chrY, chrM
- **Coordinate system:** 1-based (GFF/GTF), 0-based (BED)
- **Assembly patches:** Updates without version changes

### Model Organism References

#### Mouse
- **GRCm39/mm39:** Current mouse reference
- **GRCm38/mm10:** Previous version

#### Other Important References
- **Drosophila:** dm6 (Release 6)
- **C. elegans:** ce11 (WS245)
- **Zebrafish:** GRCz11
- **Arabidopsis:** TAIR10
- **E. coli:** NC_000913 (K-12 MG1655)

### Annotation Standards

#### Gene Annotation
- **GENCODE:** Comprehensive gene annotation (human/mouse)
- **RefSeq:** NCBI curated annotations
- **Ensembl:** Automated gene builds

#### Functional Annotation
- **GO:** Gene Ontology terms
- **KEGG:** Pathway annotations
- **Pfam:** Protein domain annotations

---

## Programming Languages & Libraries

### Python

#### Core Bioinformatics
- **Biopython:** Fundamental bioinformatics library
- **pysam:** SAM/BAM file manipulation
- **pyVCF:** VCF file parsing
- **scikit-bio:** Bioinformatics data structures

#### Data Analysis
- **pandas:** Data manipulation and analysis
- **numpy:** Numerical computing
- **scipy:** Scientific computing
- **matplotlib/seaborn:** Data visualization
- **plotly:** Interactive visualizations

#### Single-cell Analysis
- **scanpy:** Single-cell analysis toolkit
- **anndata:** Annotated data structures
- **scrublet:** Doublet detection
- **palantir:** Trajectory analysis

#### Machine Learning
- **scikit-learn:** General machine learning
- **tensorflow/pytorch:** Deep learning
- **keras:** High-level neural networks

### R

#### Core Bioconductor
- **Biostrings:** DNA/RNA/protein sequences
- **GenomicRanges:** Genomic intervals
- **IRanges:** Infrastructure for ranges
- **rtracklayer:** Import/export genomic data
- **BSgenome:** Genome sequences

#### RNA-seq Analysis
- **DESeq2:** Differential expression analysis
- **edgeR:** Count data analysis
- **limma:** Linear models for genomics
- **tximport:** Import transcript quantifications

#### Single-cell Analysis
- **Seurat:** Single-cell toolkit
- **SingleCellExperiment:** Data containers
- **scran:** Quality control and normalization
- **monocle3:** Trajectory analysis

#### Visualization
- **ggplot2:** Grammar of graphics
- **ComplexHeatmap:** Advanced heatmaps
- **circlize:** Circular visualizations
- **ggtree:** Phylogenetic tree visualization

### Other Languages

#### Perl
- **BioPerl:** Bioinformatics modules
- **Use cases:** Text processing, legacy pipelines

#### Java
- **Picard:** SAM/BAM utilities
- **GATK:** Variant calling toolkit
- **IJ (ImageJ):** Image analysis

#### C/C++
- **samtools/bcftools:** Core utilities
- **BWA/STAR:** High-performance aligners

---

## Quality Control Standards

### Sequencing Data QC

#### Read Quality Metrics
- **Per-base quality score:** Q30+ for 80%+ of bases
- **Per-sequence quality:** Mean quality > 20
- **GC content:** Should match expected distribution
- **Duplication rate:** <20% for most applications
- **Adapter contamination:** <5%

#### RNA-seq Specific
- **Alignment rate:** >85% for well-annotated genomes
- **Gene body coverage:** Uniform distribution
- **Ribosomal RNA contamination:** <10%
- **Strand specificity:** >80% for stranded protocols

### Variant Calling QC

#### Quality Filters
- **QUAL score:** >30 for high-confidence variants
- **Depth (DP):** 10-100x for WGS, 20-200x for WES
- **Genotype Quality (GQ):** >20
- **Strand bias:** Balanced read support

#### Population Filters
- **Hardy-Weinberg Equilibrium:** p > 1e-6
- **Minor Allele Frequency:** >1% for common variants
- **Call rate:** >95% across samples

### Single-cell QC

#### Cell-level Metrics
- **Number of genes:** 500-8000 per cell
- **Number of UMIs:** 1000-50000 per cell
- **Mitochondrial percentage:** <20%
- **Ribosomal percentage:** Context-dependent

#### Gene-level Metrics
- **Minimum cells:** Expressed in >3 cells
- **Maximum cells:** Not expressed in >95% of cells

---

## Common Commands & One-liners

### File Format Conversion

```bash
# SAM to BAM conversion and sorting
samtools view -b input.sam | samtools sort -o output.bam

# BAM to FASTQ extraction
samtools fastq -1 R1.fastq -2 R2.fastq input.bam

# VCF compression and indexing
bgzip variants.vcf
tabix -p vcf variants.vcf.gz

# FASTA indexing
samtools faidx reference.fasta
```

### Quality Control

```bash
# FastQC for quality assessment
fastqc *.fastq.gz -o qc_output/

# MultiQC to aggregate reports
multiqc qc_output/ -o multiqc_report/

# Basic sequence statistics
seqkit stats *.fastq.gz
```

### Alignment Commands

```bash
# BWA alignment pipeline
bwa index reference.fasta
bwa mem reference.fasta R1.fastq.gz R2.fastq.gz | \
  samtools sort -o aligned.bam
samtools index aligned.bam

# STAR alignment for RNA-seq
STAR --runMode genomeGenerate --genomeDir star_index \
     --genomeFastaFiles reference.fasta \
     --sjdbGTFfile annotation.gtf
STAR --genomeDir star_index --readFilesIn R1.fastq.gz R2.fastq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
```

### Variant Calling

```bash
# GATK variant calling
gatk HaplotypeCaller -R reference.fasta -I input.bam -O variants.vcf

# FreeBayes variant calling
freebayes -f reference.fasta input.bam > variants.vcf

# BCFtools variant calling
bcftools mpileup -f reference.fasta input.bam | \
  bcftools call -mv -o variants.vcf
```

### Text Processing

```bash
# Extract specific columns from VCF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' variants.vcf.gz

# Count reads in FASTQ
echo $(zcat file.fastq.gz | wc -l)/4 | bc

# Extract sequences by ID from FASTA
seqkit grep -p "pattern" sequences.fasta

# BED file operations
bedtools intersect -a file1.bed -b file2.bed
bedtools merge -i sorted.bed
```

### AWK One-liners

```bash
# Calculate N50 from assembly statistics
awk '{sum+=$2; sizes[NR]=$2} END {
  asort(sizes); 
  target=sum/2; 
  cumsum=0; 
  for(i=NR; i>=1; i--) {
    cumsum+=sizes[i]; 
    if(cumsum>=target) {print sizes[i]; break}
  }
}' assembly_stats.txt

# Convert GTF to BED format
awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$10}' annotation.gtf
```

---

## File Extensions Quick Reference

### Sequence Data
```
.fasta, .fa, .fas    - FASTA sequences
.fastq, .fq          - FASTQ reads with quality
.sra                 - NCBI Sequence Read Archive
.csfasta, .qual      - SOLiD colorspace format
```

### Alignment Data
```
.sam                 - Sequence Alignment/Map
.bam                 - Binary SAM
.cram                - Compressed reference alignment
.bai, .csi           - BAM/CRAM index files
```

### Variant Data
```
.vcf                 - Variant Call Format
.bcf                 - Binary VCF
.gvcf                - Genomic VCF
.tbi, .csi           - Tabix index files
```

### Annotation Data
```
.gff, .gff3          - General Feature Format
.gtf                 - Gene Transfer Format
.bed                 - Browser Extensible Data
.wig, .bigwig        - Wiggle track format
.bedgraph            - BedGraph format
```

### Analysis Results
```
.h5, .hdf5           - HDF5 hierarchical data
.h5ad                - AnnData (single-cell)
.rds                 - R data serialization
.tsv, .txt           - Tab-separated values
.csv                 - Comma-separated values
```

### Archives & Compression
```
.gz                  - Gzip compression
.bz2                 - Bzip2 compression
.tar.gz, .tgz        - Compressed archives
.zip                 - ZIP archives
```

---

## Hardware & Resource Requirements

### Memory Requirements by Analysis Type

#### Alignment
- **BWA mem:** 4-8 GB (human genome)
- **STAR:** 32-64 GB (includes genome loading)
- **HISAT2:** 8-16 GB

#### Variant Calling
- **GATK HaplotypeCaller:** 8-16 GB
- **FreeBayes:** 4-8 GB
- **Joint genotyping:** 32-64 GB (large cohorts)

#### Assembly
- **SPAdes:** 128-512 GB (mammalian genomes)
- **Canu:** 64-256 GB
- **Long-read assembly:** 512-1024 GB

#### Single-cell Analysis
- **10K cells:** 16-32 GB
- **100K cells:** 64-128 GB
- **1M+ cells:** 256-512 GB

### Storage Considerations

#### Raw Data Sizes
- **Human WGS:** 100-200 GB (FASTQ)
- **Human WES:** 10-20 GB (FASTQ)
- **RNA-seq:** 5-50 GB per sample
- **Single-cell RNA-seq:** 1-10 GB per sample

#### Processed Data
- **BAM files:** 50-100 GB (human WGS)
- **VCF files:** 1-10 GB (human genome)
- **Compressed archives:** 10-50% of original size

---

## Best Practices & Standards

### Data Management

#### File Organization
```
project/
├── data/
│   ├── raw/           # Original, unmodified data
│   ├── processed/     # Quality-controlled data
│   └── results/       # Analysis outputs
├── scripts/           # Analysis code
├── docs/              # Documentation
└── environment/       # Conda/Docker files
```

#### Naming Conventions
- Use consistent, descriptive names
- Include sample IDs, conditions, dates
- Avoid spaces and special characters
- Use snake_case or kebab-case

#### Version Control
- Use Git for code and small files
- Use Git LFS for large binary files
- Tag major analysis versions
- Document analysis parameters

### Reproducibility

#### Environment Management
```bash
# Conda environment
conda env export > environment.yml
conda env create -f environment.yml

# Docker containers
docker build -t analysis:latest .
docker run -v $(pwd):/data analysis:latest

# Package versions
pip freeze > requirements.txt
sessionInfo() # in R
```

#### Documentation Standards
- README files for each project
- Analysis notebooks with explanations
- Parameter files for reproducibility
- Data provenance tracking

### Quality Assurance

#### Analysis Checkpoints
1. **Raw data QC:** FastQC, MultiQC
2. **Alignment QC:** Mapping rates, coverage
3. **Variant QC:** Ti/Tv ratios, depth distribution
4. **Result validation:** Known positive controls

#### Statistical Considerations
- **Multiple testing correction:** FDR, Bonferroni
- **Effect size reporting:** Not just p-values
- **Confidence intervals:** For all estimates
- **Power calculations:** Sample size justification

---

*This reference guide serves as a comprehensive context document for bioinformatics workflows. It should be regularly updated as new tools, standards, and best practices emerge in the field.*
