#!/bin/bash

# Claude Code for Bioinformatics - Codespaces Setup Script
# This script sets up a complete interactive bioinformatics learning environment

set -e

echo "🚀 Setting up Claude Code for Bioinformatics Interactive Environment"
echo "=============================================================="

# Update system packages
echo "📦 Updating system packages..."
sudo apt-get update -qq

# Install essential bioinformatics tools
echo "🧬 Installing bioinformatics tools..."
sudo apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    build-essential \
    default-jre \
    python3-pip \
    python3-venv

# Install Miniconda for bioinformatics package management
echo "🐍 Installing Miniconda..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
rm miniconda.sh

# Add conda to PATH
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
export PATH="$HOME/miniconda/bin:$PATH"

# Initialize conda
$HOME/miniconda/bin/conda init bash

# Configure bioconda channels
echo "📋 Configuring bioinformatics software channels..."
$HOME/miniconda/bin/conda config --add channels defaults
$HOME/miniconda/bin/conda config --add channels bioconda
$HOME/miniconda/bin/conda config --add channels conda-forge

# Install core bioinformatics tools
echo "🔬 Installing bioinformatics software..."
$HOME/miniconda/bin/conda install -y \
    fastqc \
    multiqc \
    bwa \
    samtools \
    bcftools \
    gatk4 \
    picard \
    trimmomatic

# Install Claude Code
echo "🤖 Installing Claude Code..."
npm install -g @anthropic-ai/claude-code

# Verify Claude Code installation
echo "✅ Verifying Claude Code installation..."
claude --version || echo "⚠️ Claude Code installation may need manual verification"

# Create workspace structure
echo "📁 Setting up workspace structure..."
mkdir -p /workspaces/claude-for-bioinformatics/interactive-workspace
cd /workspaces/claude-for-bioinformatics

# Download sample data to workspace
echo "📊 Preparing sample data..."
mkdir -p interactive-workspace/sample-data
cd interactive-workspace/sample-data

# Download the educational sample data
curl -L -o sample_R1.fastq.gz "https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R1.fastq.gz"
curl -L -o sample_R2.fastq.gz "https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R2.fastq.gz"

# Validate downloads
echo "🔍 Validating sample data..."
gunzip -t sample_R1.fastq.gz && echo "✅ R1 file integrity OK" || echo "❌ R1 file may be corrupted"
gunzip -t sample_R2.fastq.gz && echo "✅ R2 file integrity OK" || echo "❌ R2 file may be corrupted"

cd /workspaces/claude-for-bioinformatics

# Set up global Claude Code context (if setup script exists)
echo "🧠 Setting up Claude Code bioinformatics context..."
if [ -f "setup.sh" ]; then
    bash setup.sh || echo "⚠️ Global setup script encountered issues"
else
    echo "ℹ️ Global setup script not found, will use manual context loading"
fi

# Create a welcome script
echo "📝 Creating welcome message..."
cat > interactive-workspace/WELCOME.md << 'EOF'
# 🎉 Welcome to Interactive Claude Code Bioinformatics Learning!

## 🚀 You're in a Real Development Environment

This is a complete bioinformatics workspace with:
- ✅ **Claude Code** installed and ready
- ✅ **Professional tools** (FastQC, MultiQC, BWA, GATK, etc.)
- ✅ **Sample data** ready for analysis
- ✅ **VS Code** with extensions for bioinformatics
- ✅ **Terminal access** for real command-line interaction

## 🔑 Next Steps: API Setup

1. **Get your Anthropic API key**:
   - Visit: https://console.anthropic.com
   - Create account or sign in
   - Go to "API Keys" and create a new key

2. **Set up your API key** (choose one method):
   
   **Method A: Environment Variable (Recommended)**
   ```bash
   export ANTHROPIC_API_KEY="your-api-key-here"
   echo 'export ANTHROPIC_API_KEY="your-api-key-here"' >> ~/.bashrc
   ```
   
   **Method B: Claude Code Config**
   ```bash
   claude auth
   # Follow prompts to enter your API key
   ```

3. **Test Claude Code**:
   ```bash
   claude --version
   claude "Hello! Can you help me with bioinformatics?"
   ```

## 📚 Interactive Tutorials Available

- **[Interactive Tutorial 1: First RNA-seq Analysis](../interactive-tutorials/01-interactive-rnaseq/README.md)**
- **[Interactive Tutorial 2: Variant Calling](../interactive-tutorials/02-interactive-variants/README.md)**
- **[Interactive Tutorial 3: Custom Workflows](../interactive-tutorials/03-custom-workflows/README.md)**

## 🎯 What Makes This Different

Unlike the Colab version where you watch code execute, here you:
- **Type actual commands** in the terminal
- **Interact with Claude Code** in real conversations
- **Build real projects** with proper file structures
- **Learn authentic workflows** used in research
- **Experience problem-solving** when things go wrong

## 🆘 Getting Help

- **Terminal issues**: Check the integrated terminal at bottom of VS Code
- **Claude Code problems**: Verify API key setup with `claude auth status`
- **Tool issues**: Check tool installation with `fastqc --version`, etc.
- **File issues**: Sample data is in `interactive-workspace/sample-data/`

**🚀 Ready to start your interactive bioinformatics journey!**
EOF

# Make executable scripts
chmod +x interactive-workspace/WELCOME.md

# Final setup message
echo ""
echo "🎉 Interactive Claude Code Bioinformatics Environment Setup Complete!"
echo ""
echo "📋 What's installed:"
echo "   ✅ Claude Code (interactive AI assistant)"
echo "   ✅ FastQC, MultiQC (quality control)"
echo "   ✅ BWA, SAMtools (alignment and processing)"
echo "   ✅ GATK4, Picard (variant calling)"
echo "   ✅ Sample RNA-seq data for learning"
echo ""
echo "🔑 Next steps:"
echo "   1. Set up your Anthropic API key (see WELCOME.md)"
echo "   2. Open interactive-workspace/WELCOME.md for full instructions"
echo "   3. Start with Interactive Tutorial 1"
echo ""
echo "🚀 Happy learning!"