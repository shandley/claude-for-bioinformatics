# Claude Code for Bioinformatics - Standard Operating Procedure

## Overview

This SOP provides a standardized approach for using Claude Code to enhance bioinformatics workflows. By following these procedures, you'll leverage AI assistance while maintaining scientific rigor and reproducibility.

**Time to implement**: 30 seconds (automated setup)  
**Immediate benefits**: Zero-friction context loading, instant project creation, consistent workflows

---

## 🚀 Quick Start Checklist (30 seconds)

### Automated Setup (Recommended)
```bash
# One command installs everything globally
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash

# Create your first project
claude-bio new rnaseq my-analysis
cd my-analysis

# Start Claude with automatic context loading
claude
```

✅ **Setup Complete** - Context documents automatically loaded, ready for analysis!

### Manual Setup (Legacy Approach)
If you prefer manual setup or the automated script doesn't work:

```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code

# Clone this repository
git clone https://github.com/shandley/claude-for-bioinformatics.git
cd claude-for-bioinformatics

# Follow manual context loading procedures below
```

---

## 📋 Standard Workflow Protocol

### Every Analysis Session Protocol

#### 1. **Context Automatically Loaded** (After Setup)
With the automated setup, bioinformatics context is automatically available!

**Automated Approach (Recommended):**
- Context documents are automatically loaded when you start Claude in any project
- No copy-pasting required
- Always up-to-date with latest versions

**Manual Approach (Legacy):**
If using manual setup, at the beginning of each Claude Code session:
```
I'm working on bioinformatics analysis. Let me provide you with domain-specific context documents to help you understand the tools, file formats, and best practices.

[Paste content from context documents - see Context Document Usage section]
```

#### 2. **Define Your Analysis Goal** (Natural Language)
Instead of trying to remember command syntax, describe what you want to accomplish:

**Good Examples**:
- "I have paired-end RNA-seq FASTQ files and need to run quality control analysis"
- "Convert these SAM files to sorted, indexed BAM files"
- "Run differential expression analysis comparing treatment vs control samples"
- "Create a variant calling pipeline following GATK best practices"

#### 3. **Review Generated Code** (Critical Safety Step)
**ALWAYS review code before execution**:
- ✅ Check file paths are correct
- ✅ Verify parameters make sense for your data
- ✅ Ensure output directories exist
- ✅ Validate that commands match your research goals

#### 4. **Execute and Validate Results**
- Run code on small test datasets first
- Check intermediate outputs
- Validate final results against known standards
- Document successful workflows for future use

#### 5. **Save Successful Workflows** (Reproducibility)
Create custom commands for repeated analyses:
```
/save-workflow "rna-seq-qc" 
# Saves the successful workflow as a reusable command
```

---

## 📚 Context Document Usage

### Essential Context Documents
Always provide these documents at session start:

#### 1. **Bioinformatics Context Reference** (Core Domain Knowledge)
```
File: bioinformatics-context-reference-guide.md
Usage: Provides Claude Code with understanding of:
- File formats (FASTQ, BAM, VCF, etc.)
- Software tools and their applications
- Quality control standards
- Best practices and conventions
```

#### 2. **Claude Code Best Practices** (Tool Usage Patterns)
```
File: claude-code-best-practices.md  
Usage: Teaches Claude Code how to:
- Structure bioinformatics projects
- Create effective custom commands
- Integrate with existing workflows
- Optimize for team collaboration
```

#### 3. **Bioinformatics One-Liners** (Command Examples)
```
File: bioinformatics-one-liners.md
Usage: Provides examples of:
- Common file manipulations
- Text processing patterns
- Quick analysis commands
- Troubleshooting approaches
```

### How to Provide Context Documents

**Option 1: Full Context (Recommended for new sessions)**
```
Here are my bioinformatics context documents:

[Paste entire content of bioinformatics-context-reference-guide.md]

[Paste entire content of claude-code-best-practices.md]

[Paste entire content of bioinformatics-one-liners.md]
```

**Option 2: Selective Context (For specific tasks)**
```
I'm working on RNA-seq analysis. Here's the relevant context:

[Paste RNA-seq sections from context documents]
```

**Option 3: Reference Context (For ongoing sessions)**
```
Please refer to the bioinformatics context I provided earlier about [specific topic]
```

---

## 🔧 Project Setup Patterns

### Standard Bioinformatics Project Structure
```
project_name/
├── data/
│   ├── raw/              # Original data files
│   ├── processed/        # Quality-controlled data
│   └── reference/        # Reference genomes, annotations
├── scripts/              # Analysis scripts and workflows
├── results/              # Analysis outputs
│   ├── qc/              # Quality control reports
│   ├── alignments/      # BAM files and indices
│   ├── variants/        # VCF files and annotations
│   └── expression/      # Gene expression results
├── docs/                 # Analysis documentation
├── .claude/
│   ├── settings.json    # Project-specific Claude Code settings
│   └── commands/        # Custom slash commands
├── CLAUDE.md            # Project context for Claude Code
└── README.md            # Project overview
```

### Essential Project Files

#### CLAUDE.md Template
Create this file in every bioinformatics project:

```markdown
# [Project Name] - Bioinformatics Analysis

## Project Overview
[Brief description of research goals and data types]

## Data Information
- **Species**: [organism]
- **Data Type**: [RNA-seq, WGS, etc.]
- **Samples**: [number and description]
- **Reference Genome**: [version used]

## Analysis Workflow
1. Quality control with FastQC/MultiQC
2. [Alignment/preprocessing steps]
3. [Primary analysis]
4. [Statistical analysis]
5. [Visualization and reporting]

## Key Commands Used
[Document successful commands for reproducibility]

## Important Notes
- [Any project-specific considerations]
- [Quality thresholds and filters applied]
- [Known issues or data quirks]
```

---

## ⚡ Common Workflow Examples

### RNA-seq Quality Control
```
Session Start: [Provide context documents]

Request: "I have paired-end RNA-seq FASTQ files in my data/raw/ directory. Please run comprehensive quality control analysis including FastQC and MultiQC reports."

Expected Output:
- FastQC analysis on all FASTQ files
- MultiQC aggregated report
- Summary of quality metrics
- Recommendations for data filtering
- Organized output in results/qc/ directory
```

### Variant Calling Pipeline
```
Session Start: [Provide context documents]

Request: "Set up a GATK best practices variant calling workflow for my human WGS samples. I have paired-end FASTQ files and want to follow current guidelines."

Expected Output:
- Complete GATK pipeline setup
- Reference genome preparation
- Alignment with BWA-MEM
- Preprocessing steps (BQSR, etc.)
- Variant calling with HaplotypeCaller
- Quality filtering recommendations
```

### File Format Conversion
```
Request: "I have various alignment files (some SAM, some unsorted BAM) and need everything converted to sorted, indexed BAM files for downstream analysis."

Expected Output:
- Automatic format detection
- Conversion pipeline (SAM→BAM)
- Coordinate sorting
- Index file generation
- Validation of output files
```

---

## 🛡️ Safety and Validation Protocols

### Before Running Any Code

#### File Safety Checklist
- [ ] **Backup important data** before analysis
- [ ] **Check available disk space** for outputs
- [ ] **Verify input file paths** are correct
- [ ] **Test on small datasets** first

#### Code Review Checklist
- [ ] **Commands make scientific sense** for your analysis
- [ ] **Parameters are appropriate** for your data type
- [ ] **Output paths are organized** and won't overwrite existing files
- [ ] **Resource requirements** (memory, CPU) are reasonable

### During Analysis

#### Monitoring Checklist
- [ ] **Check intermediate outputs** for expected file sizes
- [ ] **Monitor system resources** (memory, disk usage)
- [ ] **Validate key checkpoints** (alignment rates, variant counts)
- [ ] **Document any warnings** or unexpected results

### After Analysis

#### Results Validation
- [ ] **Compare to expected ranges** (use context document standards)
- [ ] **Check output file integrity** (not corrupted or truncated)
- [ ] **Validate against positive controls** if available
- [ ] **Document successful parameters** for future use

---

## 🔄 Team Collaboration Workflows

### Shared Project Setup
1. **Repository Structure**: Use consistent directory layout across projects
2. **Shared CLAUDE.md**: Include team-specific standards and protocols
3. **Custom Commands**: Share successful workflows as slash commands
4. **Documentation**: Maintain analysis logs and parameter records

### Code Review Process
1. **Before Implementation**: Review Claude Code outputs with team lead
2. **Parameter Validation**: Verify analysis parameters meet project standards
3. **Results Review**: Cross-check key findings with independent analysis
4. **Documentation**: Record validated workflows in shared knowledge base

### Knowledge Transfer
1. **Successful Workflows**: Document as custom slash commands
2. **Parameter Sets**: Save validated parameters for different data types
3. **Troubleshooting**: Maintain log of common issues and solutions
4. **Training**: New team members follow SOP with mentor guidance

---

## 🚨 Troubleshooting Guide

### Common Issues and Solutions

#### "Command not found" errors
**Cause**: Tool not installed or not in PATH
**Solution**: 
```
Ask Claude Code: "The analysis failed because [tool] wasn't found. Please check for alternative tools or provide installation instructions."
```

#### Memory/resource errors
**Cause**: Insufficient computational resources
**Solution**:
```
"The analysis failed due to memory limitations. Please modify the workflow to use less memory or process data in smaller batches."
```

#### File format issues
**Cause**: Unexpected file formats or corruption
**Solution**:
```
"I'm getting file format errors. Please analyze these files and recommend preprocessing steps or format conversions."
```

#### Analysis parameter questions
**Cause**: Uncertainty about appropriate parameters
**Solution**:
```
"Based on the context documents provided, what are the recommended parameters for [specific analysis] with [data type]?"
```

### Getting Help

#### From Claude Code
- Describe the error message and context
- Ask for alternative approaches
- Request explanation of analysis steps
- Seek parameter optimization advice

#### From Community
- GitHub Issues: Technical problems with workflows
- Discussion Forums: Best practices questions
- Literature: Validate approaches against published methods
- Core Facilities: Institution-specific guidance

---

## 📊 Quality Metrics and Standards

### RNA-seq Quality Thresholds
- **Per-base quality**: Q30+ for 80%+ of bases
- **Alignment rate**: >85% for well-annotated genomes
- **Duplicate rate**: <20% for most applications
- **Gene body coverage**: Uniform distribution

### Variant Calling Quality Filters
- **QUAL score**: >30 for high-confidence variants
- **Depth (DP)**: 10-100x for WGS, 20-200x for WES
- **Genotype Quality (GQ)**: >20
- **Hardy-Weinberg Equilibrium**: p > 1e-6

### Single-cell Quality Metrics
- **Genes per cell**: 500-8000 typically
- **UMIs per cell**: 1000-50000 typically
- **Mitochondrial %**: <20%
- **Doublet rate**: <10%

*Reference context documents for complete quality standards*

---

## 📈 Success Metrics

### Immediate Benefits (First Week)
- ✅ Faster analysis setup (50% time reduction)
- ✅ Fewer command-line errors
- ✅ More consistent workflows
- ✅ Better documentation

### Medium-term Benefits (First Month)
- ✅ Standardized team workflows
- ✅ Improved reproducibility
- ✅ Enhanced troubleshooting
- ✅ Knowledge transfer efficiency

### Long-term Benefits (3+ Months)
- ✅ Advanced workflow automation
- ✅ Custom tool integration
- ✅ Publication-ready analyses
- ✅ Training efficiency for new team members

---

## 🔄 Continuous Improvement

### Regular Updates
- **Monthly**: Review and update context documents
- **Quarterly**: Evaluate new tools and methods
- **As needed**: Incorporate community feedback
- **Version releases**: Update for new Claude Code features

### Feedback Collection
- **Success stories**: Document effective workflows
- **Pain points**: Identify areas for improvement
- **Team input**: Regular SOP review meetings
- **Community contributions**: Share improvements

### Documentation Maintenance
- **Analysis logs**: Maintain records of successful analyses
- **Parameter optimization**: Document best parameters for different data types
- **Tool updates**: Track version changes and compatibility
- **Training materials**: Update for new team members

---

## 📞 Support and Resources

### Internal Resources
- **Project CLAUDE.md**: Project-specific context and protocols
- **Team Documentation**: Shared workflows and standards
- **Analysis Logs**: Historical records of successful analyses
- **Mentorship**: Pair new users with experienced team members

### External Resources
- **GitHub Repository**: https://github.com/shandley/claude-for-bioinformatics
- **Context Documents**: Updated reference materials
- **Community Forum**: Questions and shared experiences
- **Issue Tracking**: Bug reports and feature requests

### Emergency Contacts
- **Primary Contact**: [Lab PI/Bioinformatics Lead]
- **Technical Support**: [Core Facility/IT Support]
- **Claude Code Support**: community forums and documentation
- **Tool-specific Support**: Official documentation and user groups

---

*This SOP is a living document. Update it based on your team's experience and evolving best practices in the field.*

**Version**: 1.0  
**Last Updated**: [Current Date]  
**Review Schedule**: Quarterly