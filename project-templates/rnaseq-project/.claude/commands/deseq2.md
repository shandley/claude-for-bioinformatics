Run differential expression analysis using DESeq2:

1. **Load count data and metadata**:
   - Import gene count matrix from results/counts/
   - Load sample metadata with experimental design information
   - Verify sample names match between counts and metadata

2. **Data preprocessing**:
   - Filter low-count genes (CPM > 1 in at least 3 samples)
   - Check for outlier samples using PCA and hierarchical clustering
   - Assess need for batch effect correction

3. **Create DESeq2 object**:
   - Set up design formula based on experimental design
   - Handle batch effects if present in metadata
   - Verify factor levels are correctly specified

4. **Run differential expression analysis**:
   - Estimate size factors for normalization
   - Estimate gene-wise dispersions
   - Fit negative binomial GLM
   - Perform Wald tests for specified contrasts

5. **Extract and filter results**:
   - Apply significance thresholds (padj < 0.05)
   - Apply fold change thresholds (|log2FC| > 1.0)
   - Generate results tables for all contrasts
   - Add gene annotations if available

6. **Quality control plots**:
   - PCA plot of samples
   - Sample-to-sample distance heatmap
   - Dispersion estimates plot
   - MA plots for each contrast
   - Volcano plots for significant genes

7. **Generate result files**:
   - Save normalized counts matrix
   - Export significant gene lists for each contrast
   - Create summary statistics table
   - Generate analysis report with plots

8. **Functional annotation preparation**:
   - Extract gene lists for GO/KEGG analysis
   - Separate up-regulated and down-regulated genes
   - Create ranked gene lists for GSEA

Use appropriate statistical parameters:
- FDR correction for multiple testing
- Cook's distance filtering for outliers
- Independent filtering for low-count genes
- Proper contrast specification for comparisons

Handle common experimental designs:
- Simple two-group comparisons
- Multi-factor designs with interaction terms
- Time course experiments
- Batch effect correction

Provide clear interpretation of results including:
- Number of significant genes per contrast
- Direction of expression changes
- Quality of the analysis (dispersion, clustering)
- Recommendations for functional analysis