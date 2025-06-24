# Module 1.1: Your First RNA-seq Analysis with Claude Code

## 🚀 Two Ways to Learn

### ⚡ **Option 1: Zero-Installation (Recommended for Beginners)**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/shandley/claude-for-bioinformatics/blob/master/guided-tutorials/01-first-rnaseq-analysis/Module_1_1_RNA_seq_Analysis.ipynb)

**Start learning immediately** - No software installation required! **[Full details →](README_Colab.md)**

### 💻 **Option 2: Local Installation (This Guide)**
**Complete tutorial below** for setting up on your own computer.

---

## Learning Objectives
By the end of this tutorial, you will:
- ✅ Successfully set up Claude Code for bioinformatics
- ✅ Create your first project using our RNA-seq template
- ✅ Run quality control analysis with guided interpretation
- ✅ Understand basic Claude Code interaction patterns
- ✅ Navigate and interpret bioinformatics output files

**Estimated Time**: 45-60 minutes  
**Prerequisites**: None (this is the starting point!)

---

## Step 1: Prerequisites Check ⚠️

### Required Reading (10 minutes)
**IMPORTANT**: Before starting this tutorial, you MUST read:
- [**Claude Code Best Practices**](../../claude-code-best-practices.md)

This covers essential setup, installation, and basic usage patterns. **You cannot skip this step.**

### Installation Verification
Ensure you have Claude Code installed:
```bash
# Check if Claude Code is installed
claude --version
```

If not installed, follow the installation instructions in the prerequisites document.

---

## Step 2: Download Sample Data (5 minutes)

For this tutorial, we'll use a small, realistic RNA-seq dataset that processes quickly while teaching real concepts.

### Get the Tutorial Dataset
```bash
# Create a workspace for this tutorial
mkdir -p ~/bioinformatics-learning/first-rnaseq
cd ~/bioinformatics-learning/first-rnaseq

# Copy sample FASTQ files from the tutorial repository
# (If you cloned the repository, the files are already available)
cp /path/to/claude-for-bioinformatics/guided-tutorials/01-first-rnaseq-analysis/sample-data/*.fastq.gz .

# OR download directly from GitHub:
curl -L -o sample_R1.fastq.gz "https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R1.fastq.gz"
curl -L -o sample_R2.fastq.gz "https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R2.fastq.gz"

# Verify files are available
ls -lh *.fastq.gz
```

**What you downloaded:**
- `sample_R1.fastq.gz`: Forward reads from a human cell line RNA-seq experiment
- `sample_R2.fastq.gz`: Reverse reads from the same experiment
- **Size**: ~50MB total (10,000 read pairs - small enough to process quickly)
- **Source**: Subset of public SRA data, quality-controlled for learning

---

## Step 3: Set Up Your Project (5 minutes)

### Run the Automated Setup
```bash
# Install bioinformatics context globally (one-time setup)
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash

# Create your first RNA-seq project
claude-bio new rnaseq first-rnaseq-tutorial
cd first-rnaseq-tutorial
```

### Understand Your Project Structure
Take a moment to explore what was created:
```bash
# Look at the project structure
tree . # or `ls -la` if tree isn't available
```

You should see:
```
first-rnaseq-tutorial/
├── data/
│   ├── raw/              # Where we'll put our FASTQ files
│   ├── processed/        # Quality-controlled data will go here
│   └── reference/        # Reference genome files
├── scripts/              # Analysis scripts
├── results/              # All outputs go here
│   ├── qc/              # Quality control reports
│   ├── alignments/      # BAM files
│   └── expression/      # Gene expression results
├── CLAUDE.md            # Project context for Claude Code
└── README.md            # Project documentation
```

### Move Your Data Into the Project
```bash
# Copy the tutorial files to the raw data directory
cp ~/bioinformatics-learning/first-rnaseq/*.fastq.gz data/raw/

# OR if you're working directly in the repository:
cp ../sample-data/*.fastq.gz data/raw/

# Verify the files are in place
ls -lh data/raw/
# Should show: sample_R1.fastq.gz (~758K) and sample_R2.fastq.gz (~811K)
```

---

## Step 4: Start Claude Code with Context (5 minutes)

### Launch Claude Code
```bash
# Start Claude Code (bioinformatics context automatically loads)
claude
```

### Verify Context Loading
When Claude Code starts, you should see confirmation that bioinformatics context was loaded. If you don't see this, make sure you ran the setup script correctly.

### Test Basic Interaction
Try this simple test to ensure everything is working:

**Type this in Claude Code:**
```
I have paired-end RNA-seq FASTQ files in my data/raw/ directory. Can you show me how to examine the basic properties of these files?
```

**Expected Response Pattern:**
Claude should provide commands to:
- Count reads in the files
- Check file compression and format
- Get basic statistics about read length and quality
- Suggest next steps for quality control

**If this works correctly**, you're ready to proceed! **If not**, check the troubleshooting section below.

---

## Step 5: Quality Control Analysis (15 minutes)

Now we'll run our first real bioinformatics analysis!

### Request Quality Control Analysis
**Type this in Claude Code:**
```
I need to run comprehensive quality control on my paired-end RNA-seq data. The files are:
- data/raw/sample_R1.fastq.gz  
- data/raw/sample_R2.fastq.gz

Please run FastQC analysis and create a MultiQC report. Organize outputs in results/qc/
```

### What Claude Should Provide
Claude should give you commands similar to:
```bash
# Create QC output directory
mkdir -p results/qc

# Run FastQC on both files
fastqc data/raw/sample_R1.fastq.gz data/raw/sample_R2.fastq.gz -o results/qc/

# Run MultiQC to combine reports
multiqc results/qc/ -o results/qc/

# Display summary
echo "Quality control complete! Check results/qc/multiqc_report.html"
```

### Execute the Commands
1. **Copy and paste** the commands Claude provided
2. **Run them one by one** in your terminal (not in Claude Code)
3. **Watch for any errors** and note what outputs are created

### Expected Processing Time
- FastQC: 1-2 minutes
- MultiQC: 10-20 seconds

---

## Step 6: Interpret Your Results (10 minutes)

### Open the QC Report
```bash
# Open the MultiQC report in your browser
open results/qc/multiqc_report.html  # macOS
# or
xdg-open results/qc/multiqc_report.html  # Linux
```

### Guided Interpretation
**Ask Claude Code to help interpret the results:**
```
I've generated FastQC and MultiQC reports. Can you help me interpret the key quality metrics? What should I be looking for in RNA-seq quality control, and what would indicate problems?
```

### Key Metrics to Understand
Claude should explain these important metrics:
1. **Per-base sequence quality**: Should be high (>28) across most of the read
2. **Per-sequence quality scores**: Most reads should have high average quality
3. **Sequence length distribution**: Should match your expected read length
4. **Duplicate levels**: Some duplication is normal in RNA-seq
5. **GC content**: Should match the expected distribution for your species
6. **Adapter content**: Should be minimal if trimming was done properly

### Hands-on Exercise
Look at your actual results and answer:
- What is the average quality score across your reads?
- Do you see any adapter contamination?
- Are there any quality issues that would require attention?
- How does your data compare to typical RNA-seq quality standards?

---

## Step 7: Next Steps Planning (5 minutes)

### Ask Claude for Workflow Guidance
```
Based on my QC results, what should be the next steps in my RNA-seq analysis workflow? I'm interested in doing differential expression analysis eventually.
```

### Expected Guidance
Claude should suggest a logical workflow:
1. **Read trimming/filtering** (if needed based on QC)
2. **Reference genome preparation**
3. **Read alignment** (STAR or HISAT2)
4. **Quantification** (featureCounts or Salmon)
5. **Differential expression** (DESeq2 or edgeR)

### Plan Your Learning Path
Based on Claude's suggestions, you can continue with:
- **Module 1.2**: Understanding Your Results (detailed interpretation)
- **Module 1.3**: Variant Calling Walkthrough (different analysis type)
- **Module 2.1**: Custom Quality Control (intermediate level)

---

## Step 8: Save Your Work (5 minutes)

### Document Your Session
Create a record of what you accomplished:

**Ask Claude:**
```
Can you help me create a summary of what we accomplished in this session? Include the commands we used and the key findings from our QC analysis.
```

### Save the Summary
Copy Claude's summary into your project documentation:
```bash
# Create a session log
echo "# First RNA-seq Analysis Session - $(date)" >> analysis_log.md
# Then paste Claude's summary into this file
```

### Commit Your Work (if using Git)
```bash
# Initialize git repository for your project
git init
git add .
git commit -m "Initial RNA-seq QC analysis - tutorial completion"
```

---

## Troubleshooting

### Common Issues and Solutions

#### "FastQC command not found"
**Problem**: FastQC is not installed
**Solution**: 
```bash
# Install FastQC (method depends on your system)
# macOS with Homebrew:
brew install fastqc

# Ubuntu/Debian:
sudo apt-get install fastqc

# Or use conda:
conda install -c bioconda fastqc
```

#### "Context not loading properly"
**Problem**: Bioinformatics context documents aren't available to Claude
**Solution**:
1. Re-run the setup script: `curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash`
2. Restart Claude Code
3. If still not working, follow the manual setup in the [Context Document Usage Guide](../../Context_Document_Usage_Guide.md)

#### "Sample data download fails"
**Problem**: Network issues or file not available
**Solution**:
1. Check your internet connection
2. Try downloading individual files
3. Use the backup sample data in this repository's sample-data directory

#### "MultiQC report is empty"
**Problem**: FastQC didn't generate expected outputs
**Solution**:
1. Check that FastQC completed successfully
2. Verify FastQC output files exist in results/qc/
3. Re-run FastQC with verbose output: `fastqc -v ...`

### Getting Help

If you encounter issues not covered here:
1. **Check the troubleshooting guide**: [Bioinformatics Troubleshooting Guide](../../context/bioinformatics-troubleshooting-guide.md)
2. **Ask Claude Code**: Describe your specific error and ask for help
3. **Community support**: Post in [GitHub Discussions](https://github.com/shandley/claude-for-bioinformatics/discussions)
4. **Report bugs**: [GitHub Issues](https://github.com/shandley/claude-for-bioinformatics/issues)

---

## What You've Accomplished

🎉 **Congratulations!** You've successfully:
- ✅ Set up Claude Code for bioinformatics research
- ✅ Created a structured analysis project
- ✅ Run your first quality control analysis
- ✅ Interpreted real bioinformatics output data
- ✅ Learned basic Claude Code interaction patterns
- ✅ Organized results in a reproducible project structure

### Skills Gained
- **Technical**: FastQC, MultiQC, file organization, command execution
- **Analytical**: Quality metric interpretation, workflow planning
- **Practical**: Claude Code usage, project management, troubleshooting

### Ready for Next Steps
You now have the foundation to:
- Continue with more advanced tutorials
- Apply these skills to your own research data
- Collaborate effectively with bioinformatics teams
- Build more complex analysis workflows

---

## Next Learning Modules

### Immediate Next Steps
- **[Module 1.2: Understanding Your Results](../02-understanding-results/)** - Deep dive into interpreting bioinformatics outputs
- **[Module 1.3: Variant Calling Walkthrough](../03-variant-calling/)** - Apply skills to a different analysis type

### When You're Ready for More
- **[Module 2.1: Custom Quality Control](../../intermediate-tutorials/01-custom-qc/)** - Customize workflows for specific needs
- **[Level 2: Independent Practice](../../intermediate-tutorials/)** - Build more complex analysis skills

---

*This tutorial is part of the Claude for Bioinformatics educational series. Your feedback helps us improve - please share your experience in our [community discussions](https://github.com/shandley/claude-for-bioinformatics/discussions)!*