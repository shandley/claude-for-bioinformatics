# Interactive Tutorial 2: Variant Calling with Claude Code

## 🎯 Learning Objectives

In this hands-on tutorial, you will:
- ✅ **Apply RNA-seq skills** to genomic variant analysis
- ✅ **Learn GATK best practices** through interactive guidance
- ✅ **Handle real genomic data** with proper quality control
- ✅ **Interpret variant calls** with biological context
- ✅ **Build variant filtering workflows** for research applications

**⏱️ Estimated Time**: 60-75 minutes  
**💰 API Cost**: ~$0.15-0.75 (minimal)  
**🔧 Prerequisites**: Complete [Interactive Tutorial 1](../01-interactive-rnaseq/) first

---

## 🚀 Before You Begin

### ✅ **Verify Your Setup**

```bash
# Check GATK and variant calling tools
gatk --version
bcftools --version
samtools --version

# Check sample data
ls ../sample-data/
```

**✋ Stop here if any tools are missing!** The Codespaces environment should have everything installed.

---

## 🧬 Step 1: Understanding Variant Calling

**💬 Start a conversation with Claude Code**:

```bash
claude
```

**💬 Begin with context**:

```
Hi Claude! I completed the RNA-seq quality control tutorial and now I want to learn variant calling. I have some genomic sequencing data and want to identify SNPs and indels using GATK best practices.

Can you explain the basic variant calling workflow and help me set up a project for this analysis?
```

**🎯 Learning focus**: Understanding the conceptual framework before diving into commands.

---

## 📊 Step 2: Data Preparation and Alignment

**💬 Continue your conversation**:

```
I have paired-end genomic sequencing data (not RNA-seq this time). The files are sample_genomic_R1.fastq.gz and sample_genomic_R2.fastq.gz in ../sample-data/.

What's the first step in preparing this data for variant calling? I've heard about BWA-MEM for alignment - can you walk me through the process?
```

**🎯 Interactive workflow**:
1. **Learn alignment basics** through conversation
2. **Set up reference genome** with Claude's guidance
3. **Run BWA-MEM alignment** step by step
4. **Understand SAM/BAM formats** through exploration

**💬 Key questions to ask**:
- "Why do we need to mark duplicates?"
- "What are read groups and why are they important?"
- "How do I know if my alignment quality is good?"

---

## 🔬 Step 3: Pre-processing for Variant Calling

**💬 Ask Claude**:

```
Now I have aligned BAM files. I understand that GATK has specific requirements for variant calling. Can you guide me through the pre-processing steps like marking duplicates and base quality score recalibration?

What's the purpose of each step and how do I know if they're working correctly?
```

**🎯 GATK workflow learning**:
1. **Mark duplicates** with Picard tools
2. **Base quality recalibration** (BQSR) with GATK
3. **Quality validation** at each step
4. **Troubleshoot issues** as they arise

**💡 Pro tip**: Ask Claude to explain the biological reasoning behind each step.

---

## 🎯 Step 4: Variant Discovery

**💬 Core variant calling conversation**:

```
Great! Now I'm ready for the actual variant calling step. I want to use GATK HaplotypeCaller to identify variants. Can you explain how this tool works and guide me through running it?

Also, what output formats should I expect and how do I interpret the initial results?
```

**🎯 Interactive discovery**:
1. **Run HaplotypeCaller** with proper parameters
2. **Understand GVCF format** and intermediate files
3. **Generate final VCF** with variant calls
4. **Initial quality assessment** of variant calls

---

## 📈 Step 5: Variant Filtering and Quality Control

**💬 Quality control conversation**:

```
I now have a VCF file with variant calls, but I know raw variant calls need filtering. Can you help me understand GATK's Variant Quality Score Recalibration (VQSR) or hard filtering approaches?

What quality metrics should I look at, and how do I decide what constitutes a high-quality variant?
```

**🎯 Quality assessment**:
1. **Learn filtering strategies** for different variant types
2. **Apply quality filters** with bcftools or GATK
3. **Generate quality metrics** and interpretation
4. **Compare before/after filtering** results

---

## 🔍 Step 6: Variant Annotation and Interpretation

**💬 Biological interpretation**:

```
Now I have filtered, high-quality variants. How do I add biological annotation to understand what these variants might mean? Can you guide me through variant annotation and help me interpret the results?

I'm particularly interested in identifying potentially functional variants.
```

**🎯 Annotation workflow**:
1. **Add functional annotations** (gene locations, effect predictions)
2. **Include population frequencies** and clinical databases
3. **Prioritize variants** by biological significance
4. **Generate summary reports** for downstream analysis

---

## 📊 Step 7: Visualization and Reporting

**💬 Results presentation**:

```
I want to create visualizations and summary reports from my variant analysis. What are the best ways to visualize variant data, and how can I create publication-quality figures?

Can you also help me generate a summary report of the key findings?
```

**🎯 Visualization skills**:
1. **Create variant density plots** and quality distributions
2. **Generate IGV sessions** for manual inspection
3. **Build summary statistics** and tables
4. **Prepare results** for sharing with collaborators

---

## 🛠️ Step 8: Automation and Workflow Development

**💬 Advanced workflow conversation**:

```
This manual process works well for learning, but I want to create an automated variant calling pipeline that I can reuse with different datasets. Can you help me build a script that implements the entire workflow?

I want something robust that includes quality checks and handles different input scenarios.
```

**🎯 Pipeline development**:
1. **Design modular workflow** with error checking
2. **Implement automated quality gates** 
3. **Add logging and progress tracking**
4. **Test pipeline** with sample data

---

## 🎓 What You've Accomplished

### 🔬 **Advanced Technical Skills**
- ✅ **Genomic data alignment** - BWA-MEM best practices
- ✅ **GATK variant calling** - Industry-standard pipeline
- ✅ **Quality control workflows** - Filtering and validation
- ✅ **Variant annotation** - Biological interpretation
- ✅ **Pipeline development** - Automated workflows

### 🧠 **Conceptual Understanding**  
- ✅ **Variant calling theory** - From reads to variants
- ✅ **Quality assessment** - Recognizing reliable variants
- ✅ **Biological interpretation** - Functional significance
- ✅ **Workflow optimization** - Efficiency and reproducibility

### 💼 **Research Applications**
- ✅ **Population genomics** - Large-scale variant analysis
- ✅ **Clinical applications** - Pathogenic variant identification
- ✅ **Comparative genomics** - Species and strain differences
- ✅ **Quality assurance** - Robust analytical practices

---

## 🚀 Next Steps

### **Continue Learning**
- **[Interactive Tutorial 3: Custom Workflows](../03-custom-workflows/)** - Build complex multi-step pipelines
- **[Advanced Genomics](../04-advanced-genomics/)** - Population analysis and structural variants
- **Apply to your research** - Use these skills with real data

### **Real Research Integration**
- **Start with pilot data** - Apply pipeline to small datasets first
- **Validate with known samples** - Use positive and negative controls
- **Collaborate with domain experts** - Discuss results with geneticists/clinicians
- **Document successful workflows** - Build institutional knowledge

---

## 💬 Troubleshooting

### **GATK Issues**
```bash
# Check GATK installation
gatk --list

# Verify reference genome format
samtools faidx reference.fa

# Check BAM file integrity
samtools quickcheck aligned.bam
```

### **Memory/Performance Issues**
```bash
# Monitor resource usage
htop

# Check available disk space
df -h

# Optimize GATK memory usage
gatk --java-options "-Xmx8g" HaplotypeCaller ...
```

### **Getting Help**
- **Ask Claude specific questions** about errors or unexpected results
- **Share exact error messages** for targeted troubleshooting
- **Check GATK documentation** for parameter details
- **Community forums** for complex technical issues

---

**🎉 Ready to master variant calling with AI assistance?**

**Start your conversation with `claude` and begin exploring genomic variants!**

*This tutorial builds real-world variant calling expertise through hands-on practice and AI guidance.*