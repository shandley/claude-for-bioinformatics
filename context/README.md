# Bioinformatics Context Documents

This directory contains domain-specific knowledge that gets automatically loaded to provide Claude Code with bioinformatics expertise.

## Context Documents

### [bioinformatics-context-reference-guide.md](bioinformatics-context-reference-guide.md)
**Purpose**: Comprehensive bioinformatics domain knowledge
- File formats (FASTQ, BAM, VCF, etc.)
- Software tools and their applications
- Quality control standards and thresholds
- Reference genomes and databases
- Best practices for different analysis types

### [bioinformatics-one-liners.md](bioinformatics-one-liners.md)
**Purpose**: Practical command examples and patterns
- Common file manipulations
- Text processing for bioinformatics data
- Tool-specific command patterns
- Troubleshooting one-liners

## How These Are Used

When you run the setup script, these documents are:
1. **Downloaded globally** to `~/.claude/bioinformatics/context/`
2. **Automatically loaded** when Claude Code starts in any project
3. **Kept up-to-date** through the update mechanism

This provides Claude Code with deep bioinformatics knowledge without requiring manual copy-pasting.

## What's NOT Here

**General Claude Code usage** is documented separately in:
- [`claude-code-best-practices.md`](../claude-code-best-practices.md) - General Claude Code usage patterns
- Project documentation and setup guides
- Non-bioinformatics specific information

This separation ensures that:
- Bioinformatics context remains focused on domain knowledge
- General Claude Code guidance can be referenced independently
- Context loading is efficient and relevant