# Claude Code for Bioinformatics

**Learn bioinformatics with AI assistance - Start analyzing data in 30 seconds!**

## 🚀 **NEW: Zero-Installation Learning**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/shandley/claude-for-bioinformatics/blob/master/guided-tutorials/01-first-rnaseq-analysis/Module_1_1_RNA_seq_Analysis.ipynb)

**👆 Click to start your first RNA-seq analysis right now - no software installation required!**

✨ **Professional bioinformatics tools (FastQC, MultiQC) in your browser**  
✨ **Real research data and publication-quality results**  
✨ **Complete tutorial with step-by-step guidance**  
✨ **Download results to continue on your computer**

[![GitHub stars](https://img.shields.io/github/stars/shandley/claude-for-bioinformatics.svg)](https://github.com/shandley/claude-for-bioinformatics/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/shandley/claude-for-bioinformatics.svg)](https://github.com/shandley/claude-for-bioinformatics/network)
[![GitHub issues](https://img.shields.io/github/issues/shandley/claude-for-bioinformatics.svg)](https://github.com/shandley/claude-for-bioinformatics/issues)

---

## 📚 Choose Your Learning Path

### 🎯 **Option 1: Instant Learning (Recommended for Beginners)**

**Start analyzing RNA-seq data in your browser immediately:**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/shandley/claude-for-bioinformatics/blob/master/guided-tutorials/01-first-rnaseq-analysis/Module_1_1_RNA_seq_Analysis.ipynb)

- ⚡ **30 seconds to start** - No downloads, no installation, no setup
- 🧬 **Real bioinformatics tools** - FastQC, MultiQC, conda environments
- 📊 **Professional results** - Generate publication-quality QC reports
- 💾 **Take results with you** - Download everything for your research
- 🎓 **Learn by doing** - See exactly how tools are installed and used

**Perfect for**: First-time users, students, anyone wanting immediate hands-on experience

### 💻 **Option 2: Local Installation (For Ongoing Research)**

**Set up Claude Code on your computer for daily research workflows:**

```bash
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash
```

**What you get:**
- **Automated global setup** - context automatically loaded in all projects
- **AI assistant with bioinformatics domain knowledge** 
- **Project creation tools** with `claude-bio new`
- **Automated context loading** - no manual document copying required
- **Quality assurance protocols** based on established field standards

**Perfect for**: Researchers ready to integrate Claude Code into daily workflows

### 🎆 **Key Features (Both Options)**
✅ Natural language interface for complex bioinformatics commands  
✅ Automated context document management  
✅ Built-in quality control standards and best practices  
✅ Structured workflow templates  
✅ Team collaboration tools  
✅ Reproducible analysis documentation  

## ⚠️ Prerequisites

### [**→ Claude Code Best Practices**](claude-code-best-practices.md) - **READ THIS FIRST**
**Essential foundation before using bioinformatics features:**
- Claude Code installation and basic setup
- Project organization and CLAUDE.md fundamentals
- Essential commands and configuration options
- Team collaboration patterns
- **⚡ Start here if you're new to Claude Code**

## 📖 Learning Resources

### 🕰️ **Learning Progression**
1. **[Claude Code Best Practices](claude-code-best-practices.md)** - Essential foundation (read first)
2. **[Zero-Installation Tutorial](guided-tutorials/01-first-rnaseq-analysis/README_Colab.md)** - Hands-on learning guide
3. **[Complete Learning Roadmap](ENHANCED_EDUCATIONAL_PLAN.md)** - Progressive skill development plan

### 🚀 **Quick Access**
- **Instant Tutorial**: Click the Colab badges above for immediate hands-on learning
- **Production Setup**: Use the installation command for ongoing research workflows
- **Advanced Learning**: Explore project templates and comprehensive guides below

## 🏆 What You'll Achieve

### 🔬 **In Your First 30 Minutes**
- ✅ **Analyze real RNA-seq data** with professional bioinformatics tools
- ✅ **Generate publication-quality reports** (FastQC and MultiQC)
- ✅ **Learn quality control interpretation** with guided explanations
- ✅ **Download results** to use in your research or presentations
- ✅ **Understand AI-assisted workflows** for future bioinformatics projects

### 📊 **Real Research Value**
- **Sample Data**: 10,000 paired-end reads with realistic quality patterns
- **Professional Tools**: Same software used in research labs worldwide
- **Publication Ready**: Generate figures and reports suitable for papers
- **Skills Transfer**: Apply techniques immediately to your own data

## 📋 Production Workflows

### [**→ Claude Code Bioinformatics SOP**](Claude_Code_Bioinformatics_SOP.md)

**Production-ready workflow guide for research teams:**
- Installation and setup checklist
- Session workflow protocols  
- Safety and validation procedures
- Team collaboration patterns
- Troubleshooting guide

## 📚 Context Documents (The Secret Sauce)

The setup script automatically provides Claude Code with bioinformatics domain expertise:

### [**→ Context Document Usage Guide**](Context_Document_Usage_Guide.md)
Complete instructions for maximizing Claude Code effectiveness (legacy approach)

### Bioinformatics Context Documents (Auto-loaded after setup)
- [**Bioinformatics Context Reference**](context/bioinformatics-context-reference-guide.md) - File formats, tools, quality standards
- [**Bioinformatics One-Liners**](context/bioinformatics-one-liners.md) - Command examples and patterns
- [**Troubleshooting Guide**](context/bioinformatics-troubleshooting-guide.md) - Error diagnosis and problem-solving
- [**Computational Resources**](context/bioinformatics-computational-resources.md) - Resource planning and optimization  
- [**Statistical Methods**](context/bioinformatics-statistical-methods.md) - Statistical approaches and best practices

### General Claude Code Documentation
- [**Claude Code Best Practices**](claude-code-best-practices.md) - General usage patterns and advanced techniques (not bioinformatics-specific)

## 🧬 Project Templates

**Copy these project structures** for instant setup:

### [RNA-seq Analysis](project-templates/rnaseq-project/)
- Complete differential expression workflow
- Quality control automation
- Custom commands included
- Publication-ready outputs

### [Variant Calling](project-templates/variant-calling-project/)  
- GATK best practices implementation
- Population genomics support
- Quality filtering protocols
- Clinical interpretation ready

### [Quality Control Analysis](project-templates/qc-analysis-project/)
- Multi-platform QC workflows
- Contamination screening
- Automated reporting
- Pass/fail assessment

## 💡 How It Works

### Before: Traditional Approach
```bash
# Remember complex syntax
bwa mem -t 8 -M -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA" \
    reference.fa read1.fq read2.fq | \
    samtools sort -@ 8 -o sample1.bam
samtools index sample1.bam
# Repeat for every sample, copy-paste context documents...
```

### After: One Setup + Natural Language
```bash
# One-time setup
curl -fsSL [setup-url] | bash

# Create project
claude-bio new rnaseq my-analysis
cd my-analysis

# Start analysis (instant context loading)
claude
> I have paired-end FASTQ files and need to align them using BWA-MEM
```
**Result**: Complete workflow with error checking, batch processing, quality validation, and automatic bioinformatics context.

## 🎯 Example Session

### Step 1: One-Time Global Setup
```bash
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash
# Downloads context documents globally, installs claude-bio helper
```

### Step 2: Create Project & Start Analysis
```bash
claude-bio new rnaseq cancer-study
cd cancer-study
claude
```

### Step 3: Describe Your Goal (Context Auto-Loaded)
```
I have RNA-seq data from a cancer study. I need to run quality control, 
align to GRCh38, and perform differential expression analysis between 
tumor and normal samples.
```

### Step 4: Get Expert Workflow (With Full Context)
Claude Code provides:
- Complete quality control pipeline
- Optimized alignment parameters
- Statistical analysis setup
- Visualization scripts
- Quality validation checks

## 🛡️ Built-in Safety

### Validation Protocols
- ✅ **Code review required** before execution
- ✅ **Quality thresholds** based on field standards  
- ✅ **Parameter validation** against best practices
- ✅ **Result verification** with known controls

### Team Standards
- 📁 **Consistent project organization**
- 📝 **Documented workflows** 
- 🔄 **Reproducible analyses**
- 👥 **Knowledge sharing** across team members

## 🌟 Key Benefits

- **Standardized workflows** ensure consistent analysis approaches across team members
- **Automated quality control** helps identify common data quality issues based on established field standards
- **Reduced learning curve** for researchers new to bioinformatics  
- **Improved reproducibility** through documented workflows and best practices
- **Natural language interface** reduces need to memorize complex command syntax

## 📖 Documentation

### For Immediate Use
- [**SOP Guide**](Claude_Code_Bioinformatics_SOP.md) - Start here for lab implementation
- [**Context Usage**](Context_Document_Usage_Guide.md) - Maximize Claude Code effectiveness
- [**Project Templates**](project-templates/) - Ready-to-use project structures

### For Future Development  
- [**Development Roadmap**](DEVELOPMENT_ROADMAP.md) - 16-week educational site development plan
- [**Content Structure**](CONTENT_STRUCTURE.md) - Complete site architecture design
- [**Educational Framework**](EDUCATIONAL_FRAMEWORK.md) - Multi-level learning progression

## 🚀 Getting Started

### Option 1: Automated Setup (Recommended)
```bash
# One command sets up everything globally
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash

# Create your first project
claude-bio new rnaseq my-analysis
cd my-analysis && claude
```

### Option 2: Manual Setup (Legacy)
1. **Install Claude Code**: `npm install -g @anthropic-ai/claude-code`
2. **Follow the [SOP Guide](Claude_Code_Bioinformatics_SOP.md)** for manual approach
3. **Use [context usage guide](Context_Document_Usage_Guide.md)** for copy-paste method

## 🤝 Contributing

We welcome contributions from the bioinformatics community!

### Ways to Contribute
- 🐛 **Report bugs** in workflows or documentation
- 💡 **Suggest improvements** to the SOP or templates
- 📝 **Share your custom commands** and project templates
- 🧪 **Test the SOP** with your data and provide feedback
- 📚 **Improve documentation** based on your experience

### How to Contribute
1. **Fork this repository**
2. **Create a feature branch** (`git checkout -b feature/amazing-workflow`)
3. **Make your changes** and test them
4. **Submit a pull request** with a clear description

## 📊 Project Status

### ✅ **Phase 1: Automated SOP** (COMPLETED)
- One-command global setup script
- Automated context document management  
- claude-bio helper command suite
- Essential project templates
- Ready for immediate lab adoption

### 🚧 **Phase 2: Educational Enhancement** (IN PROGRESS)  
- Progressive learning modules building on SOP foundation
- Hands-on tutorials with real datasets
- Guided walkthroughs for common workflows
- Community learning features

### 📅 **Phase 3: Advanced Learning** (PLANNED)
- Advanced customization tutorials
- Team collaboration workflows
- Community contribution systems
- Video demonstrations and interactive content

## 💬 Support & Community

### Getting Help
- 🐛 **Issues**: [GitHub Issues](https://github.com/shandley/claude-for-bioinformatics/issues) for bugs and feature requests
- 💬 **Discussions**: [GitHub Discussions](https://github.com/shandley/claude-for-bioinformatics/discussions) for questions and sharing
- 📧 **Direct Contact**: shandley@wustl.edu for urgent issues

### Stay Updated
- ⭐ **Star this repository** to get updates
- 👀 **Watch releases** for new features
- 🍴 **Fork and customize** for your lab's needs

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Bioinformatics Community** - For feedback and real-world testing
- **Anthropic** - For creating Claude Code and supporting computational research
- **Contributors** - Everyone who has improved this resource

## 📈 Citation

If you use this SOP in your research, please cite:

```bibtex
@software{handley2024claude_bioinformatics,
  title={Claude Code for Bioinformatics: Standard Operating Procedures and Best Practices},
  author={Handley, Scott and Contributors},
  year={2024},
  url={https://github.com/shandley/claude-for-bioinformatics},
}
```

---

**Ready to enhance your bioinformatics workflows?** 

**Quick Start**: `curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash`

**Manual Setup**: [**SOP Guide →**](Claude_Code_Bioinformatics_SOP.md)

*Built with ❤️ by the computational biology community*