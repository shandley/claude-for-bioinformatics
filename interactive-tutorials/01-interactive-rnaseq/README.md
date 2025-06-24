# Interactive Tutorial 1: Your First Claude Code RNA-seq Analysis

## 🎯 Learning Objectives

In this hands-on tutorial, you will:
- ✅ **Actually use Claude Code** in real conversations
- ✅ **Type real commands** and see authentic responses  
- ✅ **Build a complete project** from scratch with Claude's help
- ✅ **Experience genuine problem-solving** when issues arise
- ✅ **Learn interactive workflows** used in professional research

**⏱️ Estimated Time**: 45-60 minutes  
**💰 API Cost**: ~$0.10-0.50 (very minimal)  
**🔧 Prerequisites**: Anthropic API key (see WELCOME.md for setup)

---

## 🚀 Before You Begin

### ✅ **Verify Your Setup**

Open a terminal in VS Code and run these commands:

```bash
# Check Claude Code is installed
claude --version

# Check your API key is working
claude "Hello! Are you ready to help with bioinformatics?"

# Check bioinformatics tools
fastqc --version
multiqc --version

# Check sample data is available
ls ../sample-data/
```

**✋ Stop here if any of these fail!** See WELCOME.md for troubleshooting.

---

## 🎭 Interactive Learning: You Drive, Claude Guides

Unlike passive tutorials, **you'll have real conversations with Claude Code** to solve bioinformatics problems. Here's how it works:

### **Your Role**: The Researcher
- Ask Claude Code for help with your analysis goals
- Request specific commands and explanations
- Troubleshoot issues as they arise
- Make decisions about next steps

### **Claude's Role**: The AI Assistant
- Provides bioinformatics expertise and commands
- Explains concepts and reasoning
- Helps interpret results
- Suggests next steps and best practices

---

## 📁 Step 1: Project Setup with Claude Code

**🎬 Action**: Start a conversation with Claude Code to set up your project.

Open the terminal and start Claude Code:

```bash
claude
```

**💬 Have this conversation** (type these messages to Claude):

```
Hi Claude! I'm a researcher learning bioinformatics and I have some RNA-seq data I'd like to analyze. Can you help me set up a proper project structure for quality control analysis?

I have paired-end FASTQ files from an RNA-seq experiment. What's the best way to organize my project directory?
```

**📝 What to expect**: Claude should suggest a directory structure and provide commands to create it.

**🎯 Your task**: Follow Claude's suggestions to create the project structure. Ask follow-up questions if anything is unclear.

**💡 Learning note**: This is how real bioinformatics collaboration works - discussing goals and getting expert guidance before diving into commands.

---

## 🧬 Step 2: Data Exploration with Claude

**💬 Continue your conversation**:

```
Great! Now I have some sample RNA-seq data files. They're in ../sample-data/ and called sample_R1.fastq.gz and sample_R2.fastq.gz. 

Can you help me explore these files to understand what we're working with? What should I look at first?
```

**🎯 Your tasks**:
1. Copy the data files to your project as Claude suggests
2. Run the exploration commands Claude provides
3. Ask Claude to explain what you're seeing
4. **Be curious!** Ask questions about anything you don't understand

**💬 Example follow-up questions**:
- "What do these quality scores mean?"
- "Why are there four lines per read?"
- "How many reads do we have total?"
- "What does the file size tell us?"

---

## 🔬 Step 3: Quality Control Analysis

**💬 Ask Claude**:

```
Now I want to run quality control analysis on these files. I've heard of FastQC and MultiQC. Can you walk me through running a complete QC analysis and explain what we're looking for?
```

**🎯 Interactive workflow**:
1. **Get commands from Claude** for running FastQC
2. **Execute the commands** in your terminal
3. **Ask Claude to interpret** any output or errors
4. **Request next steps** based on the results
5. **Troubleshoot together** if anything goes wrong

**💡 Pro tip**: If something fails, tell Claude exactly what happened. Real bioinformatics involves lots of troubleshooting!

**💬 Example problem-solving conversation**:
```
Claude, I ran the FastQC command but got this error: [paste actual error message]

What does this mean and how can I fix it?
```

---

## 📊 Step 4: Results Interpretation

**💬 After QC analysis completes**:

```
Great! I've generated the FastQC and MultiQC reports. Can you help me understand how to interpret these results? What should I be looking for in terms of data quality?

Also, based on these results, what would be the next steps in a typical RNA-seq analysis workflow?
```

**🎯 Interactive interpretation**:
1. **Ask Claude** to explain specific metrics you see
2. **Request guidance** on whether your data quality is acceptable
3. **Discuss next steps** for downstream analysis
4. **Learn decision-making** criteria for data processing

---

## 🛠️ Step 5: Custom Workflow Development

**💬 Advanced conversation**:

```
This manual process works, but I'd like to create a reusable workflow for future RNA-seq projects. Can you help me create a script that automates the QC process and generates a summary report?

I want something I can use repeatedly with different datasets.
```

**🎯 Workflow automation**:
1. **Collaborate with Claude** to design a reusable script
2. **Build the script incrementally** with Claude's help
3. **Test the script** on your sample data
4. **Refine and improve** based on results

---

## 🎓 Step 6: Reflection and Next Steps

**💬 Wrap-up conversation**:

```
Thanks for helping me through this analysis! Can you summarize what we accomplished and suggest how I might apply these skills to my own research data?

What are the key concepts I should remember, and what would be good next steps for learning more advanced bioinformatics?
```

**🎯 Learning consolidation**:
1. **Review key concepts** with Claude
2. **Discuss real-world applications** 
3. **Plan your learning path** for advanced topics
4. **Save successful workflows** for future use

---

## 🏆 What You've Accomplished

By the end of this tutorial, you will have:

### 🔬 **Technical Skills**
- ✅ **Real Claude Code interaction** - Authentic AI-assisted workflow
- ✅ **Project organization** - Professional directory structures
- ✅ **Quality control analysis** - FastQC and MultiQC interpretation
- ✅ **Problem-solving** - Troubleshooting real bioinformatics issues
- ✅ **Workflow development** - Building reusable analysis scripts

### 🧠 **Conceptual Understanding**
- ✅ **AI-assisted research** - How to collaborate effectively with AI
- ✅ **Bioinformatics workflows** - Standard practices and decision points
- ✅ **Data quality assessment** - Critical evaluation skills
- ✅ **Professional practices** - Real-world bioinformatics approaches

### 💼 **Research Skills**
- ✅ **Interactive problem-solving** - Working through challenges collaboratively
- ✅ **Documentation habits** - Recording successful workflows
- ✅ **Continuous learning** - Asking questions and seeking explanations
- ✅ **Tool integration** - Combining AI assistance with traditional tools

---

## 🚀 Next Steps

### **Immediate Practice**
- **Try with your own data** - Apply these skills to real research
- **Experiment with variations** - Ask Claude about different approaches
- **Build your workflow library** - Save successful command patterns

### **Advanced Learning**
- **[Interactive Tutorial 2: Variant Calling](../02-interactive-variants/)** - Apply skills to genomics
- **[Interactive Tutorial 3: Custom Workflows](../03-custom-workflows/)** - Build complex pipelines
- **[Local Installation Guide](../../Claude_Code_Bioinformatics_SOP.md)** - Set up on your computer

### **Real Research Application**
- **Start small** - Use Claude Code for simple tasks in your research
- **Build confidence** - Practice with non-critical analyses first
- **Share knowledge** - Teach colleagues about AI-assisted workflows
- **Contribute back** - Share successful patterns with the community

---

## 💬 Troubleshooting

### **Claude Code Issues**
```bash
# Check API key
claude auth status

# Reset authentication
claude auth login

# Test connection
claude "test message"
```

### **Tool Issues**
```bash
# Verify tool installation
fastqc --version
multiqc --version

# Check conda environment
conda list | grep fastqc
```

### **Data Issues**
```bash
# Check sample data
ls -la ../sample-data/
file ../sample-data/*.gz
gunzip -t ../sample-data/*.gz
```

### **Getting Help**
- **Ask Claude directly** - Describe your specific issue
- **Check error messages** - Copy exact error text to Claude
- **GitHub Discussions** - Community support for complex issues
- **Documentation** - Reference guides for specific tools

---

**🎉 Ready to start your interactive bioinformatics journey?**

**Open your terminal, type `claude`, and begin the conversation!**

*Remember: This is authentic learning - you're developing real skills that you'll use in your research career.*