#!/bin/bash

# Claude Code for Bioinformatics - Automated Setup Script
# This script sets up bioinformatics context documents and helper commands globally

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
REPO_URL="https://github.com/shandley/claude-for-bioinformatics"
CONTEXT_DIR="$HOME/.claude/bioinformatics"
CLAUDE_CONFIG_DIR="$HOME/.claude"
CLAUDE_CONFIG_FILE="$CLAUDE_CONFIG_DIR/settings.json"
BIN_DIR="/usr/local/bin"
VERSION="1.0.0"

# Helper functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if running on supported system
check_system() {
    log_info "Checking system compatibility..."
    
    case "$(uname -s)" in
        Darwin*)    OS="macOS" ;;
        Linux*)     OS="Linux" ;;
        CYGWIN*|MINGW*|MSYS*) OS="Windows" ;;
        *)          OS="Unknown" ;;
    esac
    
    if [[ "$OS" == "Unknown" ]]; then
        log_error "Unsupported operating system: $(uname -s)"
        exit 1
    fi
    
    log_success "System: $OS"
}

# Check if Claude Code is installed
check_claude_code() {
    log_info "Checking Claude Code installation..."
    
    if ! command -v claude &> /dev/null; then
        log_error "Claude Code is not installed or not in PATH"
        log_info "Please install Claude Code first:"
        log_info "  npm install -g @anthropic-ai/claude-code"
        exit 1
    fi
    
    CLAUDE_VERSION=$(claude --version 2>/dev/null || echo "unknown")
    log_success "Claude Code found (version: $CLAUDE_VERSION)"
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    # Check for curl or wget
    if command -v curl &> /dev/null; then
        DOWNLOADER="curl"
    elif command -v wget &> /dev/null; then
        DOWNLOADER="wget"
    else
        log_error "Neither curl nor wget found. Please install one of them."
        exit 1
    fi
    
    # Check for git (optional, but preferred)
    if ! command -v git &> /dev/null; then
        log_warning "Git not found. Will use direct downloads instead of git clone."
        USE_GIT=false
    else
        USE_GIT=true
    fi
    
    log_success "Dependencies checked"
}

# Create directory structure
setup_directories() {
    log_info "Setting up directory structure..."
    
    # Create Claude config directory if it doesn't exist
    mkdir -p "$CLAUDE_CONFIG_DIR"
    
    # Create bioinformatics context directory
    if [[ -d "$CONTEXT_DIR" ]]; then
        log_warning "Bioinformatics context directory already exists. Backing up..."
        mv "$CONTEXT_DIR" "${CONTEXT_DIR}.backup.$(date +%Y%m%d_%H%M%S)"
    fi
    
    mkdir -p "$CONTEXT_DIR"
    mkdir -p "$CONTEXT_DIR/context"
    mkdir -p "$CONTEXT_DIR/templates"
    
    log_success "Directory structure created"
}

# Download context documents
download_context() {
    log_info "Downloading bioinformatics context documents..."
    
    cd "$CONTEXT_DIR"
    
    if [[ "$USE_GIT" == true ]]; then
        # Use git clone for better tracking
        git clone --depth 1 "$REPO_URL.git" temp_repo
        
        # Copy bioinformatics-specific context documents
        cp temp_repo/context/bioinformatics-context-reference-guide.md context/
        cp temp_repo/context/bioinformatics-one-liners.md context/
        
        # Copy project templates
        cp -r temp_repo/project-templates/* templates/
        
        # Save version info
        echo "$VERSION" > version.txt
        echo "$(date)" > last_updated.txt
        
        # Cleanup
        rm -rf temp_repo
    else
        # Use direct downloads
        BASE_URL="https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master"
        
        if [[ "$DOWNLOADER" == "curl" ]]; then
            curl -fsSL "$BASE_URL/context/bioinformatics-context-reference-guide.md" -o context/bioinformatics-context-reference-guide.md
            curl -fsSL "$BASE_URL/context/bioinformatics-one-liners.md" -o context/bioinformatics-one-liners.md
        else
            wget -q "$BASE_URL/context/bioinformatics-context-reference-guide.md" -O context/bioinformatics-context-reference-guide.md
            wget -q "$BASE_URL/context/bioinformatics-one-liners.md" -O context/bioinformatics-one-liners.md
        fi
        
        # Save version info
        echo "$VERSION" > version.txt
        echo "$(date)" > last_updated.txt
    fi
    
    log_success "Context documents downloaded"
}

# Configure Claude Code
configure_claude() {
    log_info "Configuring Claude Code for bioinformatics..."
    
    # Create or update Claude Code settings
    if [[ -f "$CLAUDE_CONFIG_FILE" ]]; then
        # Backup existing config
        cp "$CLAUDE_CONFIG_FILE" "${CLAUDE_CONFIG_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
        log_info "Existing Claude config backed up"
    fi
    
    # Create basic configuration with bioinformatics support
    cat > "$CLAUDE_CONFIG_FILE" << EOF
{
  "bioinformatics": {
    "enabled": true,
    "contextPath": "$CONTEXT_DIR/context",
    "templatesPath": "$CONTEXT_DIR/templates",
    "autoLoad": true,
    "version": "$VERSION",
    "lastUpdated": "$(date)"
  },
  "autoAcceptEdits": false,
  "theme": "dark"
}
EOF
    
    log_success "Claude Code configured"
}

# Install helper commands
install_commands() {
    log_info "Installing claude-bio helper commands..."
    
    # Check if we can write to /usr/local/bin
    if [[ ! -w "$BIN_DIR" ]] && [[ "$OS" != "Windows" ]]; then
        log_warning "Cannot write to $BIN_DIR. Helper commands will be installed to ~/.local/bin"
        BIN_DIR="$HOME/.local/bin"
        mkdir -p "$BIN_DIR"
        
        # Add to PATH if not already there
        if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
            echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
            echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc 2>/dev/null || true
            log_info "Added $BIN_DIR to PATH. Please restart your shell or run: source ~/.bashrc"
        fi
    fi
    
    # Create claude-bio command
    cat > "$BIN_DIR/claude-bio" << 'EOF'
#!/bin/bash

# Claude-Bio Helper Commands
# Manages bioinformatics projects with Claude Code

CONTEXT_DIR="$HOME/.claude/bioinformatics"
VERSION="1.0.0"

show_help() {
    cat << HELP
Claude-Bio - Bioinformatics helper for Claude Code

Usage: claude-bio <command> [options]

Commands:
  new <type> [name]     Create new bioinformatics project
                        Types: rnaseq, variants, qc
  
  list                  List available project templates
  
  update                Update bioinformatics context documents
  
  status                Check installation status and version
  
  templates             Show available project templates
  
  version               Show version information
  
  help                  Show this help message

Examples:
  claude-bio new rnaseq my-cancer-study
  claude-bio new variants population-analysis
  claude-bio update
  claude-bio status

For more information, visit:
https://github.com/shandley/claude-for-bioinformatics
HELP
}

create_project() {
    local project_type="$1"
    local project_name="${2:-${project_type}-analysis}"
    
    if [[ ! -d "$CONTEXT_DIR/templates" ]]; then
        echo "Error: Bioinformatics context not found. Please run the setup script first."
        exit 1
    fi
    
    case "$project_type" in
        rnaseq|rna-seq)
            template_dir="$CONTEXT_DIR/templates/rnaseq-project"
            ;;
        variants|variant-calling)
            template_dir="$CONTEXT_DIR/templates/variant-calling-project"
            ;;
        qc|quality-control)
            template_dir="$CONTEXT_DIR/templates/qc-analysis-project"
            ;;
        *)
            echo "Error: Unknown project type '$project_type'"
            echo "Available types: rnaseq, variants, qc"
            exit 1
            ;;
    esac
    
    if [[ ! -d "$template_dir" ]]; then
        echo "Error: Template not found: $template_dir"
        exit 1
    fi
    
    if [[ -d "$project_name" ]]; then
        echo "Error: Directory '$project_name' already exists"
        exit 1
    fi
    
    echo "Creating $project_type project: $project_name"
    cp -r "$template_dir" "$project_name"
    
    # Update CLAUDE.md to reference global context
    if [[ -f "$project_name/CLAUDE.md" ]]; then
        cat > "$project_name/CLAUDE.md.new" << NEWCLAUDEMD
# Bioinformatics Context (Auto-loaded)
The bioinformatics context documents are automatically loaded from your global Claude configuration.
This includes:
- File formats and tools reference
- Quality control standards  
- Best practices and workflows
- Command-line examples

$(cat "$project_name/CLAUDE.md")
NEWCLAUDEMD
        mv "$project_name/CLAUDE.md.new" "$project_name/CLAUDE.md"
    fi
    
    echo "Project created successfully!"
    echo "Next steps:"
    echo "  cd $project_name"
    echo "  claude"
    echo ""
    echo "The bioinformatics context will be automatically loaded."
}

check_status() {
    echo "Claude-Bio Status:"
    echo "=================="
    echo "Version: $VERSION"
    
    if [[ -d "$CONTEXT_DIR" ]]; then
        echo "Context directory: $CONTEXT_DIR ✓"
        if [[ -f "$CONTEXT_DIR/version.txt" ]]; then
            echo "Context version: $(cat "$CONTEXT_DIR/version.txt")"
        fi
        if [[ -f "$CONTEXT_DIR/last_updated.txt" ]]; then
            echo "Last updated: $(cat "$CONTEXT_DIR/last_updated.txt")"
        fi
    else
        echo "Context directory: Not found ✗"
        echo "Please run the setup script first."
    fi
    
    if command -v claude &> /dev/null; then
        echo "Claude Code: Available ✓"
    else
        echo "Claude Code: Not found ✗"
    fi
}

list_templates() {
    echo "Available Project Templates:"
    echo "============================"
    echo "rnaseq        - RNA sequencing analysis pipeline"
    echo "variants      - Variant calling with GATK best practices"  
    echo "qc            - Quality control analysis workflows"
    echo ""
    echo "Usage: claude-bio new <template> [project-name]"
}

# Main command processing
case "${1:-help}" in
    new)
        if [[ -z "$2" ]]; then
            echo "Error: Project type required"
            echo "Usage: claude-bio new <type> [name]"
            echo "Types: rnaseq, variants, qc"
            exit 1
        fi
        create_project "$2" "$3"
        ;;
    list|templates)
        list_templates
        ;;
    status)
        check_status
        ;;
    update)
        echo "Updating bioinformatics context..."
        echo "Feature coming soon. For now, re-run the setup script."
        ;;
    version)
        echo "claude-bio version $VERSION"
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        echo "Unknown command: $1"
        echo "Run 'claude-bio help' for usage information"
        exit 1
        ;;
esac
EOF
    
    chmod +x "$BIN_DIR/claude-bio"
    log_success "claude-bio command installed to $BIN_DIR"
}

# Verify installation
verify_installation() {
    log_info "Verifying installation..."
    
    # Check bioinformatics-specific context files
    local context_files=(
        "$CONTEXT_DIR/context/bioinformatics-context-reference-guide.md"
        "$CONTEXT_DIR/context/bioinformatics-one-liners.md"
    )
    
    for file in "${context_files[@]}"; do
        if [[ ! -f "$file" ]]; then
            log_error "Missing context file: $file"
            exit 1
        fi
    done
    
    # Check claude-bio command
    if command -v claude-bio &> /dev/null; then
        log_success "claude-bio command available"
    else
        log_warning "claude-bio command not in PATH. You may need to restart your shell."
    fi
    
    log_success "Installation verified"
}

# Main installation function
main() {
    echo -e "${BLUE}Claude Code for Bioinformatics Setup${NC}"
    echo "====================================="
    echo ""
    
    check_system
    check_claude_code
    check_dependencies
    setup_directories
    download_context
    configure_claude
    install_commands
    verify_installation
    
    echo ""
    echo -e "${GREEN}🎉 Installation completed successfully!${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Restart your shell or run: source ~/.bashrc"
    echo "2. Create your first project: claude-bio new rnaseq my-analysis"
    echo "3. Start Claude Code: cd my-analysis && claude"
    echo ""
    echo "The bioinformatics context will be automatically loaded!"
    echo ""
    echo "For help: claude-bio help"
    echo "For status: claude-bio status"
    echo ""
    echo "Documentation: https://github.com/shandley/claude-for-bioinformatics"
}

# Run main function
main "$@"