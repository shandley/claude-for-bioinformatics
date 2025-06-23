# Bioinformatics Statistical Methods Guide

## Table of Contents
1. [Multiple Testing Correction](#multiple-testing-correction)
2. [Differential Expression Analysis](#differential-expression-analysis)
3. [Variant Filtering and Quality Control](#variant-filtering-and-quality-control)
4. [Power Analysis and Sample Size](#power-analysis-and-sample-size)
5. [Population Genetics Statistics](#population-genetics-statistics)
6. [Single-Cell Analysis Statistics](#single-cell-analysis-statistics)
7. [Common Statistical Pitfalls](#common-statistical-pitfalls)
8. [Experimental Design Considerations](#experimental-design-considerations)

---

## Multiple Testing Correction

### The Multiple Testing Problem

When testing thousands of hypotheses simultaneously (genes, variants, etc.), the probability of false positives increases dramatically without correction.

#### Example: Why Correction is Needed
```
Testing 20,000 genes at α = 0.05
Expected false positives without correction: 20,000 × 0.05 = 1,000 genes
This is clearly unacceptable for biological interpretation
```

### False Discovery Rate (FDR) Control

#### Benjamini-Hochberg (BH) Procedure
**Most commonly used in bioinformatics**

```r
# R implementation
p_values <- c(0.001, 0.01, 0.03, 0.05, 0.08)
p_adjusted <- p.adjust(p_values, method = "BH")

# Manual calculation:
# 1. Sort p-values in ascending order
# 2. For rank i out of m tests: p_adj = p_i × m/i
# 3. Ensure monotonicity (p_adj[i] >= p_adj[i-1])
```

**When to use BH:**
- RNA-seq differential expression
- ChIP-seq peak calling
- Microarray analysis
- Any high-throughput hypothesis testing

#### Storey q-value
**More powerful than BH when many true positives exist**

```r
# Using qvalue package
library(qvalue)
q_values <- qvalue(p_values)$qvalues

# Estimates proportion of true nulls (π₀)
# More powerful when π₀ < 1 (many real effects)
```

**When to use q-value:**
- Large-scale genomics studies
- When expecting many true positives
- Meta-analyses with strong signals

### Family-Wise Error Rate (FWER) Control

#### Bonferroni Correction
**Most conservative approach**

```r
p_bonferroni <- p.adjust(p_values, method = "bonferroni")
# p_adj = min(1, p × m)
```

**When to use Bonferroni:**
- Small number of tests (< 100)
- When any false positive is unacceptable
- Candidate gene studies
- Clinical diagnostic tests

#### Holm-Bonferroni
**Less conservative than Bonferroni**

```r
p_holm <- p.adjust(p_values, method = "holm")
```

### Practical Guidelines

#### Choosing Correction Method
```
Number of Tests | Recommended Method
< 10           | No correction needed
10-100         | Bonferroni or Holm
100-1,000      | BH-FDR
> 1,000        | BH-FDR or q-value
```

#### Implementation Examples

**DESeq2 (RNA-seq):**
```r
# Uses BH-FDR by default
results(dds, alpha = 0.05)  # padj column is BH-adjusted
```

**GATK Variant Filtering:**
```bash
# No multiple testing correction typically applied
# Quality filters used instead (QUAL, QD, etc.)
```

**ChIP-seq Peak Calling:**
```r
# MACS2 uses q-value (FDR)
macs2 callpeak -q 0.05  # q-value threshold
```

---

## Differential Expression Analysis

### Count-Based Methods (RNA-seq)

#### DESeq2
**Recommended for most RNA-seq analyses**

```r
library(DESeq2)

# Model: Negative binomial with empirical Bayes shrinkage
dds <- DESeq(dds)  # Estimates size factors, dispersions, fits model

# Key assumptions:
# 1. Count data follows negative binomial distribution
# 2. Few genes are differentially expressed
# 3. Library sizes can be normalized

# When to use:
# - Standard RNA-seq (bulk)
# - Small to medium sample sizes (n < 100 per group)
# - When biological replicates are available
```

#### edgeR
**Alternative to DESeq2, similar performance**

```r
library(edgeR)

# Uses empirical Bayes to moderate dispersions
y <- DGEList(counts)
y <- calcNormFactors(y)  # TMM normalization
design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)

# When to use:
# - Similar to DESeq2
# - Some prefer for time-course experiments
# - When you want more control over normalization
```

#### limma-voom
**Good for complex experimental designs**

```r
library(limma)

# Transforms count data for linear modeling
v <- voom(counts, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# When to use:
# - Complex experimental designs
# - Time-course experiments
# - When you're familiar with limma
```

### Microarray Analysis

#### limma
**Standard for microarray differential expression**

```r
# For continuous expression values
fit <- lmFit(expression_matrix, design)
fit <- eBayes(fit)
topTable(fit, adjust="BH")
```

### Effect Size Considerations

#### Fold Change Thresholds
```r
# Common thresholds:
# |log2FC| > 1.0  (2-fold change) - Standard
# |log2FC| > 0.58 (1.5-fold change) - Liberal
# |log2FC| > 1.58 (3-fold change) - Conservative

# DESeq2 with fold change threshold
results(dds, alpha=0.05, lfcThreshold=1.0)
```

#### Log Fold Change Shrinkage
```r
# DESeq2 lfcShrink for more reliable effect sizes
res_shrink <- lfcShrink(dds, coef="condition_treated_vs_control", type="ashr")

# Benefits:
# - Reduces noise in log fold changes
# - Better for ranking genes
# - More reliable for downstream analysis
```

### Sample Size Considerations

#### Minimum Replicates
```
Analysis Type    | Minimum Reps | Recommended
Pilot study      | 2-3         | 3-4
Standard DE      | 3           | 5-6
Complex design   | 4-5         | 6-8
Clinical study   | 6-8         | 10-12
```

#### Power Analysis for RNA-seq
```r
# Using RNASeqPower package
library(RNASeqPower)

# Estimate sample size needed
rnapower(depth=10, n=5, cv=0.4, effect=1.5, alpha=0.05)

# Parameters:
# depth: average read depth per gene
# cv: coefficient of variation
# effect: fold change to detect
# alpha: significance level
```

---

## Variant Filtering and Quality Control

### Single Variant Quality Metrics

#### GATK Quality Scores
```bash
# Key metrics and typical thresholds:

# QUAL: Phred-scaled quality score
# Threshold: QUAL > 30 (99.9% confidence)

# QD: QualByDepth (QUAL/DP)
# SNPs: QD > 2.0
# Indels: QD > 2.0

# FS: FisherStrand (strand bias)
# SNPs: FS < 60.0
# Indels: FS < 200.0

# MQ: RMSMappingQuality
# Threshold: MQ > 40.0

# SOR: StrandOddsRatio
# SNPs: SOR < 3.0
# Indels: SOR < 10.0
```

#### Depth-Based Filtering
```bash
# Depth (DP) filtering considerations:
# Too low: unreliable genotype calls
# Too high: repetitive regions, copy number variants

# Typical thresholds:
# WGS: 10 ≤ DP ≤ 100
# WES: 20 ≤ DP ≤ 200

# Population-specific considerations:
# Adjust upper bound based on coverage distribution
```

### Population-Level Quality Control

#### Hardy-Weinberg Equilibrium
```bash
# Test for HWE departure
# p-value threshold: p > 10⁻⁶

# Violations may indicate:
# - Genotyping errors
# - Population stratification
# - Inbreeding
# - Copy number variants
# - Natural selection

# Implementation:
vcftools --vcf input.vcf --hardy --out hwe_stats
awk '$8 < 1e-6 {print}' hwe_stats.hwe  # Variants violating HWE
```

#### Allele Frequency Filtering
```bash
# Minor Allele Frequency (MAF) filtering:
# Removes rare variants that may be errors

# Common thresholds:
# MAF > 0.01 (1%) - Standard
# MAF > 0.05 (5%) - Conservative
# MAF > 0.001 (0.1%) - Liberal

vcftools --vcf input.vcf --maf 0.01 --recode --out filtered
```

#### Missing Data Filtering
```bash
# Individual call rate
# Remove samples with >10% missing genotypes
vcftools --vcf input.vcf --missing-indv --out missing_stats

# Variant call rate  
# Remove variants missing in >5% of samples
vcftools --vcf input.vcf --max-missing 0.95 --recode --out filtered
```

### Variant Recalibration (VQSR)

#### When to Use VQSR
```
Sample Size | Recommendation
< 30 exomes | Hard filtering
≥ 30 exomes | VQSR
< 10 genomes| Hard filtering  
≥ 10 genomes| VQSR
```

#### VQSR Training Resources
```bash
# High-confidence training sets:
# HapMap 3.3 sites (SNPs) - sensitivity: 99.7%
# Omni 2.5M chip (SNPs) - sensitivity: 99.7%
# 1000G Phase 3 (SNPs) - sensitivity: 99.0%
# Mills and 1000G (indels) - sensitivity: 99.0%

# Truth sensitivity thresholds:
# SNPs: 99.7% (retains ~99.7% of true variants)
# Indels: 99.0% (more conservative due to higher error rate)
```

---

## Power Analysis and Sample Size

### Principles of Power Analysis

#### Key Components
```
Power = 1 - β (Type II error rate)
α = Type I error rate (significance level)
Effect size = magnitude of difference to detect
Sample size = number of observations per group

Standard values:
Power: 0.80 (80%)
α: 0.05 (5%)
```

### RNA-seq Power Analysis

#### Effect Size Considerations
```r
# Biological vs Statistical significance
# |log2FC| > 1.0: 2-fold change (often biologically meaningful)
# |log2FC| > 0.58: 1.5-fold change (may be meaningful)

# Factors affecting power:
# 1. Read depth per gene
# 2. Biological variability (CV)
# 3. Number of replicates
# 4. Effect size to detect
```

#### Sample Size Calculation
```r
library(RNASeqPower)

# For 80% power to detect 2-fold change
rnapower(depth=20, cv=0.4, effect=2, alpha=0.05, power=0.8)

# Typical results:
# High expression genes: n=3-4 per group
# Medium expression genes: n=5-6 per group  
# Low expression genes: n=8-12 per group
```

### GWAS Power Analysis

#### Sample Size Requirements
```r
# For common variants (MAF > 5%):
# Small effect (OR = 1.2): n > 10,000 cases
# Medium effect (OR = 1.5): n > 1,000 cases
# Large effect (OR = 2.0): n > 200 cases

# For rare variants (MAF < 1%):
# Much larger sample sizes needed
# Consider gene-based tests
```

#### Tools for Power Calculation
```r
# QUANTO: epidemiological studies
# CATS: case-control studies  
# GAS Power Calculator: online tool
# PAWE-3D: three-dimensional power analysis
```

### Single-Cell Power Analysis

#### Unique Considerations
```r
# Power depends on:
# 1. Number of cells per condition
# 2. Number of UMIs per cell
# 3. Percentage of cells expressing gene
# 4. Effect size (log fold change)

# Typical requirements:
# Well-expressed genes: 100-500 cells per condition
# Lowly-expressed genes: 1000+ cells per condition
```

---

## Population Genetics Statistics

### Measures of Genetic Diversity

#### Nucleotide Diversity (π)
```
π = average number of nucleotide differences per site between sequences
Typical values:
- Humans: π ≈ 0.001 (1 difference per 1000 bp)
- Drosophila: π ≈ 0.01
- Plants: highly variable (0.001-0.02)
```

#### Watterson's θ
```
θ = 4Neμ (for diploids)
Ne = effective population size
μ = mutation rate per site per generation

Estimator: θw = S/Σ(1/i) where S = number of segregating sites
```

#### Tajima's D
```
D = (π - θw) / √Var(π - θw)

Interpretation:
D = 0: neutral evolution
D < 0: population expansion or purifying selection
D > 0: population bottleneck or balancing selection

Typical range: -2 to +2
```

### Population Structure

#### FST (Fixation Index)
```
FST = (HT - HS) / HT
HT = total heterozygosity
HS = subpopulation heterozygosity

Interpretation:
0.00-0.05: little differentiation
0.05-0.15: moderate differentiation  
0.15-0.25: great differentiation
>0.25: very great differentiation
```

#### Principal Component Analysis (PCA)
```r
# Using SNPRelate
library(SNPRelate)

# Convert VCF to GDS format
snpgdsVCF2GDS(vcf.fn="input.vcf", out.fn="output.gds")

# Perform PCA
pca <- snpgdsPCA(genofile, autosome.only=FALSE)

# Plot results
plot(pca$eigenvect[,1], pca$eigenvect[,2], 
     xlab="PC1", ylab="PC2")

# Variance explained:
pc.percent <- pca$varprop * 100
```

### Natural Selection Tests

#### McDonald-Kreitman Test
```
Tests for adaptive evolution by comparing:
- Synonymous vs non-synonymous substitutions (between species)
- Synonymous vs non-synonymous polymorphisms (within species)

Neutrality index (NI) = (Pn/Ps) / (Dn/Ds)
NI < 1: adaptive evolution
NI > 1: purifying selection
```

#### Selective Sweep Detection
```r
# Extended Haplotype Homozygosity (EHH)
# Integrated Haplotype Score (iHS)
# Cross-population Extended Haplotype Homozygosity (XP-EHH)

# Using rehh package
library(rehh)
ihs_results <- scan_hh(haplohh_data, limhaplo=2, limehh=0.05)
```

---

## Single-Cell Analysis Statistics

### Quality Control Metrics

#### Cell-Level QC
```r
# Key metrics and typical thresholds:

# Number of genes per cell
# Low: < 500 (poor quality cells)
# High: > 8000 (potential doublets)
# Typical range: 1000-6000

# Number of UMIs per cell  
# Low: < 1000 (poor quality)
# High: > 50000 (potential doublets)
# Typical range: 2000-30000

# Mitochondrial gene percentage
# High: > 20% (stressed/dying cells)
# Typical: < 10-15%

# Implementation in Seurat:
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
```

#### Gene-Level QC
```r
# Genes expressed in minimum number of cells
# Typical threshold: genes in ≥ 3 cells
# This removes noise while retaining rare cell type markers

# Implementation:
counts <- GetAssayData(object = seurat_obj, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 3
```

### Normalization Methods

#### Log Normalization
```r
# Standard approach: log(CPM/100 + 1)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

# Assumptions:
# - Captures relative differences well
# - Simple and interpretable
# - May not handle zeros optimally
```

#### SCTransform
```r
# Newer method addressing UMI sampling noise
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt")

# Benefits:
# - Models UMI sampling noise explicitly  
# - Better handling of highly variable genes
# - Integrated variance stabilization
```

#### Size Factor Normalization
```r
# Using scran package
library(scran)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)

# Benefits:
# - Robust to composition effects
# - Good for differential expression
```

### Dimensionality Reduction

#### Principal Component Analysis
```r
# Select highly variable genes first
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                  nfeatures = 2000)

# Scale data and run PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Determine significant PCs
ElbowPlot(seurat_obj, ndims = 50)
# Typically use 10-50 PCs for downstream analysis
```

#### UMAP and t-SNE
```r
# Non-linear dimensionality reduction for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# Parameters to consider:
# - Number of PCs to include
# - Perplexity (t-SNE): typical range 10-100
# - n.neighbors (UMAP): typical range 5-50
```

### Clustering and Cell Type Identification

#### Graph-Based Clustering
```r
# Build k-nearest neighbor graph
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Find clusters using Louvain algorithm
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Resolution parameter:
# Lower (0.1-0.5): fewer, larger clusters
# Higher (0.8-2.0): more, smaller clusters
# Optimize based on biological expectations
```

#### Marker Gene Identification
```r
# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                         min.pct = 0.25, logfc.threshold = 0.25)

# Statistical test options:
# "wilcox": Wilcoxon rank sum (default, non-parametric)
# "bimod": Likelihood-ratio test
# "roc": ROC analysis
# "t": Student's t-test
# "negbinom": Negative binomial
```

### Differential Expression in Single Cells

#### Considerations Unique to scRNA-seq
```r
# Challenges:
# 1. High dropout rates (excess zeros)
# 2. Overdispersion
# 3. Batch effects
# 4. Cell type composition effects

# Methods addressing these:
# - MAST: Hurdle model for zero-inflation
# - scDD: Handles different types of differential expression
# - ZINB-WaVE: Zero-inflated negative binomial
```

#### Pseudobulk Analysis
```r
# Aggregate cells by sample, then use bulk RNA-seq methods
# Often more powerful than single-cell specific methods

library(muscat)
pb <- aggregateData(sce, assay = "counts", fun = "sum",
                   by = c("cluster_id", "sample_id"))

# Then use DESeq2 or edgeR on pseudobulk data
```

---

## Common Statistical Pitfalls

### Multiple Comparisons Issues

#### Hidden Multiple Testing
```r
# Common scenarios where multiple testing is overlooked:

# 1. Subgroup analyses
# Testing treatment effect in multiple patient subgroups
# Solution: Plan subgroup analyses, adjust for multiple comparisons

# 2. Multiple time points
# Testing differential expression at each time point separately
# Solution: Use time-course specific methods

# 3. Multiple conditions
# Pairwise comparisons between all conditions
# Solution: Use omnibus test first, then post-hoc comparisons
```

#### Cherry-Picking Results
```r
# Problems:
# - Reporting only significant results
# - Post-hoc hypotheses based on data
# - Multiple analysis approaches until significant

# Solutions:
# - Pre-specify hypotheses and analysis plan
# - Report all tested hypotheses
# - Use proper multiple testing correction
```

### Sample Size and Power Issues

#### Post-Hoc Power Analysis
```r
# Problem: Calculating power after seeing results
# "Observed power" is misleading - directly related to p-value

# Better approach:
# - Calculate power for effect sizes of biological interest
# - Use confidence intervals to assess precision
# - Consider effect size estimation uncertainty
```

#### Unbalanced Designs
```r
# Problems with severely unbalanced groups:
# - Reduced power
# - Violates assumptions of some tests
# - Difficult interpretation

# Solutions:
# - Plan for balanced designs when possible
# - Use methods robust to imbalance
# - Consider matching or stratification
```

### Normalization and Batch Effects

#### Inappropriate Normalization
```r
# RNA-seq: Don't use RPKM/FPKM for differential expression
# Use TMM, DESeq2 size factors, or similar methods

# Single-cell: Consider composition effects
# Cells with very different expression profiles can bias normalization

# Variants: Don't normalize quality scores
# Use raw quality metrics for filtering decisions
```

#### Batch Effect Confounding
```r
# Problem: Batch effects confounded with biological variables
# Example: All cases in batch 1, all controls in batch 2

# Solutions:
# - Randomize samples across batches
# - Include batch as covariate if possible
# - Use batch correction methods (Combat, Limma removeBatchEffect)
# - Report batch information in publications
```

### Correlation vs Causation

#### Gene Expression Correlations
```r
# High correlation doesn't imply:
# - Direct regulation
# - Causal relationship
# - Functional interaction

# Better approaches:
# - Time-course data for temporal relationships
# - Perturbation experiments
# - Mendelian randomization for causal inference
# - Network analysis with multiple data types
```

### Statistical Model Assumptions

#### Normality Assumptions
```r
# Many tests assume normal distributions
# RNA-seq: Count data is not normal (use negative binomial models)
# Solution: Use appropriate models for data type

# Check assumptions:
qqnorm(residuals); qqline(residuals)
shapiro.test(residuals)  # For small samples
```

#### Independence Assumptions
```r
# Violations common in:
# - Repeated measures (multiple time points per individual)
# - Family studies (related individuals)
# - Spatial data (neighboring samples)

# Solutions:
# - Mixed effects models
# - Family-based tests
# - Account for correlation structure
```

---

## Experimental Design Considerations

### Randomization and Controls

#### Proper Randomization
```r
# Randomize at appropriate level:
# - Treatment assignment
# - Sample processing order
# - Sequencing lane assignment
# - Analysis order

# Block randomization for known confounders:
# - Age groups
# - Sex
# - Batch processing dates
```

#### Control Types
```r
# Negative controls:
# - Untreated samples
# - Vehicle-only treatment
# - Mock transfection

# Positive controls:
# - Known responsive genes/pathways
# - Spike-in controls for sequencing
# - Technical replicates

# Internal controls:
# - Housekeeping genes (though these can vary)
# - Spike-in RNA for normalization
```

### Replication Strategies

#### Biological vs Technical Replicates
```r
# Biological replicates:
# - Different individuals/samples
# - Measures biological variability
# - Required for statistical inference

# Technical replicates:
# - Same sample, processed multiple times
# - Measures technical variability
# - Less important with modern protocols

# Recommendation: Prioritize biological replicates
```

#### Minimum Sample Sizes
```r
# General guidelines:

# RNA-seq differential expression:
# Pilot: n=3 per group (minimum)
# Standard: n=5-6 per group
# Robust: n=8-10 per group

# GWAS:
# Common variants: n=1000+ cases (minimum)
# Rare variants: n=10,000+ cases

# Single-cell:
# Cell type identification: 100-1000 cells per type
# Differential expression: 500+ cells per condition
```

### Time-Course Experiments

#### Design Considerations
```r
# Number of time points:
# Minimum: 3-4 time points
# Better: 6+ time points for smooth trajectories

# Time point selection:
# - Include baseline (time 0)
# - Expected peak response time
# - Return to baseline (if applicable)
# - Early and late response phases

# Statistical approaches:
# - Spline-based methods
# - Time-course specific DE methods (maSigPro, STEM)
# - Trajectory analysis (monocle, slingshot)
```

### Multi-Factor Designs

#### Factorial Designs
```r
# 2x2 example: Treatment × Genotype
# Main effects: Treatment effect, Genotype effect
# Interaction: Treatment effect differs by genotype

# Design matrix in R:
design <- model.matrix(~ treatment * genotype)

# Benefits:
# - Test multiple hypotheses efficiently
# - Detect interactions
# - More generalizable results

# Considerations:
# - Sample size requirements increase
# - Interpretation becomes more complex
```

#### Confounding Variables
```r
# Common confounders in genomics:
# - Age, sex, ancestry
# - Batch effects
# - Cell type composition
# - Environmental factors

# Strategies:
# - Matching on confounders
# - Including as covariates
# - Stratified analysis
# - Propensity score methods
```

---

*This statistical methods guide provides frameworks for making informed decisions about analysis approaches in bioinformatics. Always consider the specific assumptions and limitations of chosen methods in the context of your experimental design and research questions.*