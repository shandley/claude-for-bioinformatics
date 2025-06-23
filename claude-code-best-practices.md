# Claude Code Best Practices - Lab SOP Guide

## Table of Contents
1. [Setup & Configuration](#setup--configuration)
2. [Essential Files & Dotfiles](#essential-files--dotfiles)
3. [Core Features & Commands](#core-features--commands)
4. [Project Structure & Organization](#project-structure--organization)
5. [Advanced Techniques](#advanced-techniques)
6. [Team Collaboration](#team-collaboration)
7. [Troubleshooting & Best Practices](#troubleshooting--best-practices)
8. [Cost Optimization](#cost-optimization)

---

## Setup & Configuration

### Initial Setup
1. **Install Claude Code**
   ```bash
   npm install -g @anthropic-ai/claude-code
   ```

2. **Authentication Setup**
   - Obtain API key from Anthropic Console
   - Run `claude` and follow authentication prompts
   - API keys are stored securely (macOS Keychain on Mac)

3. **Essential Configuration Commands**
   ```bash
   # Terminal setup (enables Shift+Enter for newlines, bell notifications)
   /terminal-setup
   
   # IDE integration (connects to VS Code, reads linter warnings)
   /ide
   
   # Initialize project documentation
   /init
   ```

### Configuration Files Locations

#### Global Configuration (`~/.claude/`)
- `settings.json` - Main configuration file
- `commands/` - Personal slash commands (available across all projects)
- `.claude.json` - Legacy config (being deprecated in favor of settings.json)

#### Project Configuration (`.claude/`)
- `settings.json` - Project-specific settings (checked into git)
- `settings.local.json` - Personal preferences (git-ignored)
- `commands/` - Project-specific slash commands
- `CLAUDE.md` - Project knowledge base and instructions

#### MCP Configuration
- `.mcp.json` - Project-level MCP servers (checked into git)
- `~/.claude/settings.json` - Global MCP server configuration

---

## Essential Files & Dotfiles

### CLAUDE.md - Project Knowledge Base
**Purpose:** Central hub for project-specific information, coding standards, and frequently used commands.

**Key Sections:**
```markdown
# Project Name

## Overview
Brief description of the project and its purpose.

## Development Commands
```bash
# Build
npm run build

# Test
npm test

# Lint
npm run lint
```

## Coding Standards
- Use TypeScript for all new code
- Follow ESLint configuration
- Write tests for all new features
- Use conventional commits

## Architecture Notes
- Key files and their purposes
- Important patterns and conventions
- Common gotchas and solutions

## Team Preferences
- Naming conventions
- Preferred libraries and tools
- Code review guidelines
```

**Advanced CLAUDE.md Features:**
- **Import other files:** Add `@path/to/file.md` to load additional documentation
- **Team collaboration:** Share consistent instructions across team members
- **Personal preferences:** Can include individual coding style preferences

### Custom Slash Commands (`.claude/commands/`)

#### Project Commands Example Structure:
```
.claude/commands/
├── optimize.md              # /project:optimize
├── fix-issue.md             # /project:fix-issue
├── frontend/
│   └── component.md         # /project:frontend:component
└── backend/
    └── api.md               # /project:backend:api
```

#### Sample Command Files:

**optimize.md:**
```markdown
Analyze the performance of this code and suggest three specific optimizations:
```

**fix-issue.md:**
```markdown
Find and fix issue #$ARGUMENTS. Follow these steps:
1. Use `gh issue view` to get the issue details
2. Understand the problem described in the issue
3. Search the codebase for relevant files
4. Implement the necessary changes to fix the issue
5. Write and run tests to verify the fix
6. Ensure code passes linting and type checking
7. Create a descriptive commit message
8. Push and create a PR

Remember to use the GitHub CLI (`gh`) for all GitHub-related tasks.
```

#### Personal Commands (`~/.claude/commands/`)
Commands available across all projects:

**security-review.md:**
```markdown
Review this code for security vulnerabilities, focusing on:
- Input validation
- Authentication and authorization
- SQL injection risks
- XSS vulnerabilities
- Sensitive data exposure
```

### Settings Configuration Examples

#### Global Settings (`~/.claude/settings.json`)
```json
{
  "theme": "dark-daltonized",
  "autoAcceptEdits": false,
  "mcpServers": {
    "filesystem": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-filesystem", "/path/to/workspace"]
    }
  }
}
```

#### Project Settings (`.claude/settings.json`)
```json
{
  "allowedTools": ["edit", "run", "git"],
  "testCommand": "npm test",
  "lintCommand": "npm run lint",
  "buildCommand": "npm run build"
}
```

---

## Core Features & Commands

### Essential Slash Commands
- `/help` - Show all available commands
- `/clear` - Clear conversation history (use frequently!)
- `/config` - Manage settings
- `/permissions` - Manage tool permissions
- `/mcp` - Manage MCP servers
- `/ide` - Connect to IDE
- `/compact` - Summarize conversation to reduce token usage
- `/approved-tools` - Manage tool permissions
- `/continue` - Resume previous conversation
- `/resume` - Pick specific conversation to resume

### File Operations
```bash
# Tab completion for files
@<tab>  # Shows file picker

# Reading files
Read the logging.py file
Read all files in the src/ directory

# File references
@src/components/Button.tsx  # Direct file reference
```

### Image & Visual Support
Claude Code supports multiple ways to work with images:

1. **Drag and drop** images directly into the terminal
2. **Copy/paste** screenshots (Mac: `cmd+ctrl+shift+4` to clipboard, then `ctrl+v`)
3. **File references** - specify image file paths
4. **Screenshots for UI development** - iterate on designs until they match mocks

### Extended Thinking Mode
Trigger deeper reasoning with these phrases:
- `think` - 4,000 tokens of thinking
- `think harder`, `think intensely`, `think longer` - Enhanced thinking
- `megathink` - 10,000 tokens of thinking
- `ultrathink` - 31,999 tokens of thinking

**Example Usage:**
```
Think deeply about the best approach for implementing this authentication system in our codebase.
```

### Git Integration
Claude Code integrates seamlessly with Git:
```bash
# Commit with standard format
git commit -m "$(cat <<'EOF'
feat: add user authentication

- Implement JWT token validation
- Add login/logout endpoints
- Include password hashing

🤖 Generated with Claude Code
Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Project Structure & Organization

### Recommended Directory Structure
```
project-root/
├── .claude/
│   ├── settings.json          # Project settings
│   ├── settings.local.json    # Personal preferences (git-ignored)
│   └── commands/              # Project-specific commands
│       ├── deploy.md
│       ├── test.md
│       └── docs.md
├── .mcp.json                  # MCP server configuration
├── CLAUDE.md                  # Main project documentation
├── .gitignore                 # Include .claude/settings.local.json
└── [your project files]
```

### Dotfiles Management
Store Claude Code configurations in your dotfiles repository:

```bash
# Add to dotfiles
~/.claude/settings.json
~/.claude/commands/
```

Link dotfiles in new environments:
```bash
ln -s ~/dotfiles/.claude/settings.json ~/.claude/settings.json
ln -s ~/dotfiles/.claude/commands ~/.claude/commands
```

---

## Advanced Techniques

### MCP Server Integration
**Model Context Protocol (MCP)** extends Claude Code with additional tools.

#### Popular MCP Servers:
- **Filesystem:** File system operations
- **Puppeteer:** Browser automation and screenshots
- **Brave Search:** Web search capabilities
- **GitHub:** Repository operations
- **Sentry:** Error monitoring integration

#### MCP Configuration Example:
```json
{
  "mcpServers": {
    "puppeteer": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-puppeteer"]
    },
    "brave-search": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-brave-search"],
      "env": {
        "BRAVE_API_KEY": "your-api-key"
      }
    }
  }
}
```

### Automation & Headless Mode
For CI/CD integration and automation:

```bash
# Headless mode with prompt
claude -p "Fix all linting errors" --output-format json

# Stream JSON output for real-time processing
claude -p "Analyze test coverage" --output-format stream-json

# Bypass permissions for automation
claude --dangerously-skip-permissions -p "Run tests and commit fixes"
```

### Advanced Workflow Patterns

#### Test-Driven Development (TDD)
1. Ask Claude to write failing tests first
2. Implement code to make tests pass
3. Refactor while keeping tests green
4. Claude excels at this pattern

#### Visual Development Workflow
1. Provide design mockups (drag & drop images)
2. Ask Claude to implement the design
3. Use Puppeteer MCP for screenshots
4. Iterate until design matches exactly

#### Subagent Pattern
For complex problems, use subagents:
```
Create a subagent to analyze the database schema and suggest optimizations.
Then create another subagent to review the security implications.
```

---

## Team Collaboration

### Shared Project Setup
1. **Check in essential files:**
   - `.claude/settings.json` (project settings)
   - `.claude/commands/` (shared commands)
   - `CLAUDE.md` (project documentation)
   - `.mcp.json` (MCP server configuration)

2. **Git ignore personal files:**
   ```gitignore
   .claude/settings.local.json
   ```

3. **Document team conventions in CLAUDE.md:**
   - Code style preferences
   - Testing requirements
   - Deployment procedures
   - Common commands and workflows

### Team Naming Conventions
Add personality to your team's CLAUDE.md files:
- Choose fun codenames for team members
- Create engaging, memorable project personas
- Example: "When you need help, ask MR BEEF" or "Harp Dog knows best"

### Shared Commands Library
Create a library of useful commands that benefit the entire team:

**Common Team Commands:**
- `/project:deploy` - Deployment procedures
- `/project:test-full` - Complete test suite
- `/project:security-scan` - Security review
- `/project:performance-check` - Performance analysis
- `/project:docs-update` - Documentation updates

---

## Troubleshooting & Best Practices

### Performance Optimization

#### Use /clear Frequently
AI agents become unpredictable with long conversations:
```bash
/clear  # Use this often, especially when switching tasks
```

#### Manage Context Window
```bash
/compact  # Summarize conversation to reduce token usage
```

#### Auto-Accept Edits
Enable auto-accept for routine tasks:
```bash
# Press Shift+Tab twice to enable auto-accept mode
# Useful for: linting fixes, formatting, routine refactoring
```

### Common Issues & Solutions

#### "Command not found" errors
- Ensure Claude Code is in your PATH
- Check global npm installation
- Verify authentication status

#### Hidden files not appearing with @ shortcut
- Known limitation: dotfiles (.env, .github) don't show in @ completion
- Use explicit file paths as workaround

#### Slow performance on large edits
- Break down large tasks into smaller chunks
- Use /clear to start fresh conversations
- Consider using subagents for complex multi-part tasks

#### Permission prompts
- Use `/permissions` to manage domain allowlists
- Consider `--dangerously-skip-permissions` for trusted automation

### Security Considerations

#### Safe Usage Patterns:
- Review all code changes before accepting
- Use auto-accept only for trusted, routine operations
- Regularly audit generated code for security issues
- Implement pre-commit hooks for additional validation

#### Sensitive Data:
- Never include API keys or secrets in CLAUDE.md
- Use environment variables for sensitive configuration
- Review commits to ensure no sensitive data was included

---

## Cost Optimization

### Token Management Strategies

#### Efficient Conversation Patterns:
1. **Start fresh frequently** - Use `/clear` to avoid token accumulation
2. **Be specific** - Precise prompts reduce back-and-forth
3. **Use /compact** - Summarize long conversations
4. **Batch similar tasks** - Complete related work in single sessions

#### Model Selection:
- **Primary model:** Claude 3.7+ or Claude 4 for complex reasoning
- **Haiku model:** Used automatically for lightweight operations
- Monitor usage through your Anthropic console

#### Cost-Effective Practices:
```bash
# Use headless mode for automation (no interactive tokens)
claude -p "Fix lint errors" --output-format json

# Batch file operations
Read all TypeScript files and suggest consistency improvements

# Use CLAUDE.md to provide context once rather than repeatedly
```

### Monitoring Usage
- Check Anthropic Console for token usage
- Set up spending alerts if available
- Average estimated cost: ~$6/day for typical usage

---

## Quick Reference

### Essential Setup Checklist
- [ ] Install Claude Code via npm
- [ ] Authenticate with Anthropic API
- [ ] Run `/terminal-setup` for better terminal integration
- [ ] Create project CLAUDE.md file with `/init`
- [ ] Set up essential slash commands in `.claude/commands/`
- [ ] Configure MCP servers if needed
- [ ] Add team conventions to CLAUDE.md
- [ ] Set up dotfiles management for configuration

### Daily Workflow Commands
```bash
claude                    # Start Claude Code in current directory
/clear                   # Clear conversation (use frequently!)
/ide                     # Connect to IDE for linter integration
@<tab>                   # File picker with tab completion
/project:<command>       # Use custom project commands
/compact                 # Reduce token usage
```

### Emergency Commands
```bash
/help                    # Show all available commands
Ctrl+C                   # Cancel current operation
Esc                      # Stop current task
/approved-tools          # Reset tool permissions
claude --version         # Check Claude Code version
```

---

*This SOP guide is a living document. Update it as new features and best practices emerge in your lab's Claude Code usage.*
