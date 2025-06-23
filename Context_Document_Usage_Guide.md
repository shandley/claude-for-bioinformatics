# Context Document Usage Guide

## How to Use These Documents with Claude Code

This guide explains how to effectively provide the bioinformatics context documents to Claude Code for optimal assistance.

**Note**: With the automated setup, this manual approach is no longer needed! This guide is maintained for legacy/manual setups.

## ⚠️ Prerequisites

### [**→ Claude Code Best Practices**](claude-code-best-practices.md) - **READ FIRST**
**Before using bioinformatics context documents, learn Claude Code fundamentals:**
- Basic installation and setup
- Project organization principles
- CLAUDE.md file creation and usage
- Essential commands and workflows
- **These basics are required for effective bioinformatics usage**

---

## 📋 **Quick Reference**

### At Every Session Start:
1. **Copy and paste all three context documents** (see templates below)
2. **Describe your analysis goal** in natural language
3. **Reference specific sections** during analysis
4. **Save successful workflows** for future use

---

## 📚 **Context Document Templates**

### **Session Start Template** (Copy-Paste Ready)

```
I'm working on bioinformatics analysis. Let me provide you with domain-specific context documents to help you understand the tools, file formats, and best practices in computational biology.

DOCUMENT 1: BIOINFORMATICS CONTEXT REFERENCE
===========================================
[Paste entire content of bioinformatics-context-reference-guide.md here]

DOCUMENT 2: BIOINFORMATICS ONE-LINERS
====================================
[Paste entire content of bioinformatics-one-liners.md here]

Now I'm ready to discuss my bioinformatics analysis needs with full context.
```

### **Targeted Context Template** (For Specific Analyses)

```
I'm working on [RNA-seq/variant calling/quality control] analysis. Here's the relevant bioinformatics context:

RELEVANT SECTIONS FROM CONTEXT DOCUMENTS:
[Paste specific sections related to your analysis type]

My specific analysis goal: [Describe what you want to accomplish]
```

---

## 🎯 **Document-Specific Usage**

### **1. bioinformatics-context-reference-guide.md**

**Primary Purpose**: Provides Claude Code with comprehensive domain knowledge

**Key Sections for Different Analyses**:

#### For RNA-seq Analysis:
```
Relevant sections to highlight:
- File Formats → FASTQ, FASTA (lines 19-32)
- RNA-seq Analysis Tools → STAR, Salmon, DESeq2 (lines 141-152)  
- Quality Control Standards → RNA-seq specific (lines 371-376)
- Programming Languages → R/Bioconductor (lines 317-341)
```

#### For Variant Calling:
```
Relevant sections to highlight:
- Alignment Formats → SAM/BAM/CRAM (lines 34-46)
- Variant Formats → VCF/BCF (lines 48-54)
- Variant Calling Tools → GATK, FreeBayes (lines 128-140)
- Quality Control → Variant calling QC (lines 377-389)
```

#### For Single-Cell Analysis:
```
Relevant sections to highlight:
- Other Important Formats → HDF5, AnnData (lines 76-84)
- Single-cell Analysis Tools → scanpy, Seurat (lines 154-167)
- Quality Control → Single-cell QC (lines 390-401)
```

**Usage Pattern**:
```
"Based on the bioinformatics context reference I provided, what file format should I use for storing [specific data type]?"

"According to the quality control standards in my context document, what thresholds should I use for [specific QC metric]?"
```

### **2. bioinformatics-one-liners.md**

**Primary Purpose**: Provides specific command examples and patterns

**Key Sections by Task**:
- File conversions: Lines 406-421
- Quality control: Lines 423-434  
- Alignment: Lines 436-451
- Variant calling: Lines 453-465
- Text processing: Lines 467-500

**Usage Pattern**:
```
"Using the one-liner examples I provided, help me convert SAM files to sorted BAM files."

"Based on the text processing patterns in my context, extract specific columns from this VCF file."
```

---

## 🔄 **Progressive Context Strategies**

### **Strategy 1: Full Context (Recommended for Complex Analyses)**
Provide all three documents at session start for comprehensive assistance.

**Best for**:
- New analysis pipelines
- Complex multi-step workflows  
- Training new team members
- Troubleshooting difficult problems

### **Strategy 2: Targeted Context (For Focused Tasks)**
Provide specific sections relevant to your immediate task.

**Best for**:
- Quick file conversions
- Specific tool usage questions
- Parameter optimization
- Debugging specific steps

### **Strategy 3: Reference Context (For Ongoing Sessions)**
Reference previously provided context during the session.

**Example**:
```
"Referring to the quality control standards I provided earlier, what should I do about samples with >30% duplication rate?"
```

---

## 💡 **Effective Context Phrases**

### **Starting a Session**:
- "I'm working on bioinformatics analysis. Here are my context documents..."
- "Let me provide domain-specific context for computational biology..."
- "Here's comprehensive bioinformatics reference material to guide our analysis..."

### **Referencing During Analysis**:
- "Based on the file format guide I provided..."
- "According to the quality standards in my context..."
- "Using the tool recommendations from my reference documents..."
- "Following the workflow patterns I shared..."

### **For Complex Decisions**:
- "Given the context documents I provided, what's the best approach for..."
- "Considering the standards and tools in my reference materials..."
- "Based on all the bioinformatics context I've shared..."

---

## 🎯 **Context Optimization Tips**

### **1. Highlight Relevant Sections**
When providing context, emphasize sections most relevant to your task:

```
I'm working on RNA-seq analysis. Here are my context documents, with RNA-seq sections particularly relevant:

[Full context documents]

Please pay special attention to:
- RNA-seq tools and workflows
- Quality control standards for RNA-seq
- File format handling for sequencing data
```

### **2. Provide Analysis-Specific Context**
Include details about your specific analysis:

```
Context for my analysis:
- Species: Human (GRCh38)
- Data: Paired-end RNA-seq, 50M reads per sample
- Goal: Differential expression between treatment/control
- Available tools: [list your installed software]

[Context documents follow...]
```

### **3. Update Context as Needed**
If analysis requirements change during the session:

```
"My analysis requirements have evolved. Now I also need to [new requirement]. Based on the context I provided, how should I modify the workflow?"
```

---

## 🛡️ **Quality Assurance with Context**

### **Validation Questions to Ask**:
- "Does this workflow follow the best practices outlined in my context documents?"
- "Are these parameters consistent with the quality standards I provided?"
- "Based on my context reference, are there any steps I'm missing?"
- "Do these file formats align with the standards in my documentation?"

### **Cross-Reference Checks**:
- Compare suggested tools against your context document recommendations
- Verify quality thresholds match provided standards
- Ensure workflow patterns follow established best practices
- Confirm file formats are appropriate for downstream analysis

---

## 📊 **Measuring Context Effectiveness**

### **Signs Your Context is Working Well**:
✅ Claude Code suggests appropriate tools for your data type  
✅ Recommended parameters align with published best practices  
✅ Workflow steps follow logical bioinformatics progression  
✅ Quality control suggestions match field standards  
✅ File format recommendations are suitable for downstream analysis  

### **Signs You Need Better Context**:
❌ Suggestions seem generic rather than bioinformatics-specific  
❌ Recommended tools don't match your analysis type  
❌ Quality thresholds seem arbitrary or inappropriate  
❌ Workflow steps are missing important bioinformatics considerations  
❌ File format suggestions don't align with field standards  

---

## 🔄 **Context Maintenance**

### **Regular Updates**:
- **Monthly**: Review context documents for outdated information
- **New Projects**: Add project-specific context to standard documents
- **Tool Updates**: Include new tools and version changes
- **Standards Evolution**: Update quality thresholds and best practices

### **Team Synchronization**:
- **Shared Context**: Ensure all team members use consistent context documents
- **Local Additions**: Add lab-specific protocols and standards
- **Version Control**: Track changes to context documents
- **Training**: Update new team member onboarding with current context

---

## 📞 **Troubleshooting Context Issues**

### **Problem**: Claude Code doesn't seem to understand bioinformatics concepts
**Solution**: Ensure you've provided the full bioinformatics-context-reference-guide.md

### **Problem**: Suggested workflows don't follow best practices
**Solution**: Emphasize the quality standards sections in your context

### **Problem**: Tool recommendations aren't appropriate
**Solution**: Highlight relevant tool sections and specify your available software

### **Problem**: File format suggestions are incorrect
**Solution**: Reference specific file format sections when asking questions

---

*This guide ensures you get maximum value from the context documents and enables Claude Code to provide expert-level bioinformatics assistance.*