#!/bin/bash

# Claude Code for Bioinformatics - Uninstall Script
# This script removes bioinformatics context documents and helper commands

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
CONTEXT_DIR="$HOME/.claude/bioinformatics"
CLAUDE_CONFIG_DIR="$HOME/.claude"
CLAUDE_CONFIG_FILE="$CLAUDE_CONFIG_DIR/settings.json"
BIN_DIRS=("/usr/local/bin" "$HOME/.local/bin")

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

# Confirm uninstallation
confirm_uninstall() {
    echo -e "${YELLOW}This will remove:${NC}"
    echo "  - Bioinformatics context documents ($CONTEXT_DIR)"
    echo "  - claude-bio helper command"
    echo "  - Bioinformatics configuration from Claude settings"
    echo ""
    echo "Your projects and Claude Code installation will remain unchanged."
    echo ""
    read -p "Are you sure you want to continue? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Uninstallation cancelled."
        exit 0
    fi
}

# Remove context directory
remove_context() {
    log_info "Removing bioinformatics context documents..."
    
    if [[ -d "$CONTEXT_DIR" ]]; then
        # Create backup before removal
        backup_dir="${CONTEXT_DIR}.removed.$(date +%Y%m%d_%H%M%S)"
        mv "$CONTEXT_DIR" "$backup_dir"
        log_success "Context directory moved to: $backup_dir"
        log_info "You can safely delete this backup if you don't need it"
    else
        log_warning "Context directory not found: $CONTEXT_DIR"
    fi
}

# Remove helper commands
remove_commands() {
    log_info "Removing claude-bio helper commands..."
    
    local removed=false
    for bin_dir in "${BIN_DIRS[@]}"; do
        local claude_bio_path="$bin_dir/claude-bio"
        if [[ -f "$claude_bio_path" ]]; then
            rm "$claude_bio_path"
            log_success "Removed: $claude_bio_path"
            removed=true
        fi
    done
    
    if [[ "$removed" == false ]]; then
        log_warning "claude-bio command not found in standard locations"
    fi
}

# Clean Claude configuration
clean_claude_config() {
    log_info "Cleaning bioinformatics configuration from Claude settings..."
    
    if [[ -f "$CLAUDE_CONFIG_FILE" ]]; then
        # Create backup
        backup_file="${CLAUDE_CONFIG_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
        cp "$CLAUDE_CONFIG_FILE" "$backup_file"
        
        # Remove bioinformatics section using python/jq if available
        if command -v python3 &> /dev/null; then
            python3 -c "
import json
import sys

try:
    with open('$CLAUDE_CONFIG_FILE', 'r') as f:
        config = json.load(f)
    
    if 'bioinformatics' in config:
        del config['bioinformatics']
        
        with open('$CLAUDE_CONFIG_FILE', 'w') as f:
            json.dump(config, f, indent=2)
        
        print('Bioinformatics configuration removed')
    else:
        print('No bioinformatics configuration found')
        
except Exception as e:
    print(f'Error updating config: {e}')
    sys.exit(1)
"
        elif command -v jq &> /dev/null; then
            jq 'del(.bioinformatics)' "$CLAUDE_CONFIG_FILE" > "${CLAUDE_CONFIG_FILE}.tmp"
            mv "${CLAUDE_CONFIG_FILE}.tmp" "$CLAUDE_CONFIG_FILE"
            log_success "Bioinformatics configuration removed using jq"
        else
            log_warning "Cannot automatically clean configuration (python3/jq not found)"
            log_info "Config backup created: $backup_file"
            log_info "You may manually remove the 'bioinformatics' section from $CLAUDE_CONFIG_FILE"
        fi
    else
        log_warning "Claude configuration file not found: $CLAUDE_CONFIG_FILE"
    fi
}

# Check for remaining bioinformatics projects
check_projects() {
    log_info "Checking for bioinformatics projects..."
    
    # Look for CLAUDE.md files that might reference bioinformatics context
    local found_projects=()
    
    # Search in common project locations
    for search_dir in "$HOME" "$HOME/Projects" "$HOME/Documents" "$HOME/research" "$(pwd)"; do
        if [[ -d "$search_dir" ]]; then
            while IFS= read -r -d '' file; do
                if grep -l "bioinformatics" "$file" &>/dev/null; then
                    found_projects+=("$(dirname "$file")")
                fi
            done < <(find "$search_dir" -maxdepth 3 -name "CLAUDE.md" -print0 2>/dev/null)
        fi
    done
    
    if [[ ${#found_projects[@]} -gt 0 ]]; then
        log_warning "Found potential bioinformatics projects:"
        printf '  %s\n' "${found_projects[@]}" | head -10
        if [[ ${#found_projects[@]} -gt 10 ]]; then
            echo "  ... and $((${#found_projects[@]} - 10)) more"
        fi
        echo ""
        log_info "These projects may need manual cleanup of context references"
    fi
}

# Main uninstallation function
main() {
    echo -e "${BLUE}Claude Code for Bioinformatics Uninstaller${NC}"
    echo "==========================================="
    echo ""
    
    confirm_uninstall
    
    echo ""
    log_info "Starting uninstallation..."
    
    remove_context
    remove_commands
    clean_claude_config
    check_projects
    
    echo ""
    echo -e "${GREEN}✅ Uninstallation completed successfully!${NC}"
    echo ""
    echo "What was removed:"
    echo "  - Bioinformatics context documents"
    echo "  - claude-bio helper command"
    echo "  - Bioinformatics configuration from Claude settings"
    echo ""
    echo "What remains unchanged:"
    echo "  - Claude Code installation"
    echo "  - Your existing projects"
    echo "  - Other Claude configurations"
    echo ""
    echo "If you reinstall later, run:"
    echo "  curl -fsSL https://raw.githubusercontent.com/shandley/claude-for-bioinformatics/master/setup.sh | bash"
    echo ""
    echo "Thank you for using Claude Code for Bioinformatics!"
}

# Run main function
main "$@"