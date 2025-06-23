# Troubleshooting Guide: First RNA-seq Analysis Tutorial

## Common Issues and Solutions

### Setup Issues

#### 1. "claude command not found"
**Symptoms**: Getting `command not found` when running `claude --version`
**Cause**: Claude Code is not installed or not in PATH
**Solutions**:
```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code

# Verify installation
claude --version

# If still not working, check your PATH
echo $PATH
```

#### 2. "curl command fails during setup"
**Symptoms**: Setup script download fails or hangs
**Cause**: Network connectivity or permissions issues
**Solutions**:
```bash
# Try manual download
wget https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh
chmod +x setup.sh
./setup.sh

# Or check your internet connection
ping github.com
```

#### 3. "claude-bio command not found after setup"
**Symptoms**: Setup completes but `claude-bio` command doesn't work
**Cause**: Command not in PATH or shell needs restart
**Solutions**:
```bash
# Restart your shell
exec $SHELL

# Or manually add to PATH
export PATH="$HOME/.local/bin:$PATH"

# Check if command exists
which claude-bio
```

### Data Download Issues

#### 4. "Sample data download fails"
**Symptoms**: Cannot download tutorial FASTQ files
**Cause**: Network issues or file permissions
**Solutions**:
```bash
# Check if URLs are accessible
curl -I https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R1.fastq.gz

# Try alternative download method
wget -O sample_R1.fastq.gz "https://github.com/shandley/claude-for-bioinformatics/raw/master/guided-tutorials/01-first-rnaseq-analysis/sample-data/sample_R1.fastq.gz"

# Create dummy files for testing (if downloads fail)
# This creates small test files to practice with
head -4000 /path/to/any/fastq/file > sample_R1.fastq
gzip sample_R1.fastq
```

#### 5. "Downloaded files are corrupted"
**Symptoms**: FASTQ files won't open or give format errors
**Cause**: Incomplete download or file corruption
**Solutions**:
```bash
# Check file integrity
gunzip -t sample_R1.fastq.gz
gunzip -t sample_R2.fastq.gz

# Check file sizes (should be ~25MB each)
ls -lh sample_*.fastq.gz

# Re-download if corrupted
rm sample_*.fastq.gz
# Run download commands again
```

### Analysis Tool Issues

#### 6. "FastQC not installed"
**Symptoms**: `fastqc: command not found`
**Cause**: FastQC software not available
**Solutions**:
```bash
# Install via package manager
# macOS with Homebrew:
brew install fastqc

# Ubuntu/Debian:
sudo apt-get update
sudo apt-get install fastqc

# CentOS/RHEL:
sudo yum install fastqc

# Using conda (all platforms):
conda install -c bioconda fastqc

# Verify installation
fastqc --version
```

#### 7. "MultiQC not installed"
**Symptoms**: `multiqc: command not found`
**Cause**: MultiQC Python package not available
**Solutions**:
```bash
# Install via pip
pip install multiqc

# Using conda
conda install -c bioconda multiqc

# Verify installation
multiqc --version
```

#### 8. "FastQC runs but produces no output"
**Symptoms**: FastQC command completes but no reports generated
**Cause**: Output directory doesn't exist or permission issues
**Solutions**:
```bash
# Create output directory
mkdir -p results/qc

# Check permissions
ls -la results/

# Run FastQC with verbose output
fastqc -v data/raw/sample_R1.fastq.gz -o results/qc/

# Check for error messages in output
```

#### 9. "MultiQC report is empty"
**Symptoms**: MultiQC runs but report shows no data
**Cause**: No FastQC reports found or wrong directory
**Solutions**:
```bash
# Check that FastQC outputs exist
ls -la results/qc/
# Should see files ending in _fastqc.html and _fastqc.zip

# Run MultiQC with verbose output
multiqc results/qc/ -o results/qc/ -v

# Check MultiQC log for specific errors
```

### Claude Code Interaction Issues

#### 10. "Claude doesn't understand bioinformatics terms"
**Symptoms**: Claude gives generic responses instead of bioinformatics-specific help
**Cause**: Bioinformatics context not properly loaded
**Solutions**:
```bash
# Check if setup script completed successfully
claude-bio status

# Re-run setup if needed
curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash

# Manually load context if automatic loading fails
# Follow instructions in Context_Document_Usage_Guide.md
```

#### 11. "Claude provides wrong file paths"
**Symptoms**: Commands reference incorrect directories or files
**Cause**: Claude needs project-specific context
**Solutions**:
- Always run Claude from your project directory
- Provide clear file paths in your requests
- Use `ls` commands to show Claude your actual file structure
```bash
# Example of providing context to Claude:
ls -la data/raw/
# Then ask: "I have these FASTQ files in data/raw/. Please run FastQC on them."
```

#### 12. "Generated commands don't work"
**Symptoms**: Commands from Claude fail when executed
**Cause**: Version differences, missing dependencies, or system-specific issues
**Solutions**:
- Check that all required tools are installed
- Verify file paths exist
- Look for system-specific command variations
- Ask Claude to explain what each command does
- Test commands on small datasets first

### File and Directory Issues

#### 13. "Permission denied errors"
**Symptoms**: Cannot write to output directories or execute scripts
**Cause**: Insufficient file permissions
**Solutions**:
```bash
# Check current permissions
ls -la

# Fix directory permissions
chmod 755 results/
chmod 755 results/qc/

# For script execution issues
chmod +x script_name.sh
```

#### 14. "No space left on device"
**Symptoms**: Operations fail with disk space errors
**Cause**: Insufficient disk space for outputs
**Solutions**:
```bash
# Check disk usage
df -h

# Clean up temporary files
rm -rf /tmp/*

# Use smaller test datasets
# Compress outputs when possible
```

### Operating System Specific Issues

#### 15. **macOS**: "Developer tools not installed"
**Symptoms**: Compilation errors or missing build tools
**Solutions**:
```bash
# Install Xcode command line tools
xcode-select --install

# Or install Homebrew first
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### 16. **Windows**: "Bash commands don't work"
**Symptoms**: Unix commands not recognized in Windows
**Solutions**:
- Use Windows Subsystem for Linux (WSL)
- Install Git Bash or similar Unix environment
- Use Windows-equivalent commands where possible

#### 17. **Linux**: "Package manager permissions"
**Symptoms**: Cannot install packages with apt/yum
**Solutions**:
```bash
# Use sudo for system package installation
sudo apt-get install fastqc

# Or use user-level package managers
conda install -c bioconda fastqc
pip install --user multiqc
```

## Getting Additional Help

### 1. Check Log Files
Many tools generate log files that contain detailed error information:
```bash
# FastQC logs
ls results/qc/*_fastqc.zip
# Extract and check summary.txt

# MultiQC logs
cat results/qc/multiqc_data/multiqc.log
```

### 2. Test with Minimal Examples
If complex analysis fails, try simpler versions:
```bash
# Test FastQC with one file
fastqc data/raw/sample_R1.fastq.gz

# Test with uncompressed data
gunzip -c data/raw/sample_R1.fastq.gz | head -1000 > test.fastq
fastqc test.fastq
```

### 3. Verify Tool Versions
Different versions may have different requirements:
```bash
fastqc --version
multiqc --version
claude --version
```

### 4. Community Resources
- **GitHub Issues**: [Report bugs](https://github.com/shandley/claude-for-bioinformatics/issues)
- **GitHub Discussions**: [Ask questions](https://github.com/shandley/claude-for-bioinformatics/discussions)
- **Bioinformatics Stack Exchange**: General bioinformatics questions
- **Tool-specific documentation**: FastQC, MultiQC official docs

### 5. Alternative Approaches
If one method doesn't work, try alternatives:
- Different quality control tools (fastp instead of FastQC)
- Alternative installation methods (conda instead of package manager)
- Cloud-based analysis platforms
- Docker containers for consistent environments

---

## Prevention Tips

### Before Starting
1. **System Check**: Verify all required tools are installed
2. **Disk Space**: Ensure adequate storage for outputs
3. **Permissions**: Check write access to working directories
4. **Network**: Verify internet connectivity for downloads

### During Analysis
1. **Test Small**: Always test with small datasets first
2. **Check Outputs**: Verify each step produces expected results
3. **Save Progress**: Document successful commands and parameters
4. **Monitor Resources**: Watch disk space and memory usage

### Documentation
1. **Keep Notes**: Record what works and what doesn't
2. **Version Info**: Document tool versions and system details
3. **Error Messages**: Save complete error messages for troubleshooting
4. **Success Recipes**: Document successful workflows for future use

---

*If your issue isn't covered here, please [open a GitHub issue](https://github.com/shandley/claude-for-bioinformatics/issues) with details about your system, the exact commands you ran, and the complete error messages.*