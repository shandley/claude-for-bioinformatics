# Claude for Bioinformatics - Educational Resource

## Project Overview
Development of a comprehensive educational GitHub Pages site that teaches bioinformaticians how to effectively leverage Claude Code for computational biology workflows, analysis automation, and team collaboration.

## Repository Information
- **GitHub Repository**: https://github.com/shandley/claude-for-bioinformatics
- **Remote Origin**: Connected and configured for development
- **GitHub Pages URL**: https://shandley.github.io/claude-for-bioinformatics (when deployed)

## Project Goals

### Primary Objectives
1. **Create the definitive educational resource** for using Claude Code in bioinformatics research
2. **Serve multiple skill levels** from novice computational biologists to experienced bioinformaticians
3. **Provide practical, copy-paste ready examples** that solve real-world research problems
4. **Establish best practices** for AI-assisted bioinformatics workflows
5. **Build a community resource** that research teams can adopt and customize

### Target Audience
- **Novice bioinformaticians** learning computational approaches
- **Experienced researchers** wanting to integrate AI assistance
- **Core facility managers** seeking workflow standardization
- **Graduate students** and postdocs in computational biology
- **Research teams** looking to improve collaboration and reproducibility

## Development Workflow

### Essential Development Commands

```bash
# Start development server (when Jekyll is set up)
bundle exec jekyll serve --drafts

# Build site for production
bundle exec jekyll build

# Git workflow
git add .
git commit -m "feat: add new content section"
git push origin main

# Preview changes locally
open http://localhost:4000
```

### Content Development Process
1. **Plan content in planning documents** (see DEVELOPMENT_ROADMAP.md)
2. **Write content in markdown** following educational framework
3. **Test examples locally** to ensure they work
4. **Review for clarity and accuracy** across skill levels
5. **Deploy to GitHub Pages** for community feedback

### Quality Standards
- **All code examples must be tested and functional**
- **Clear explanations for every concept introduced**
- **Progressive disclosure** - simple concepts before complex ones
- **Real-world context** for every example and workflow
- **Accessible language** avoiding unnecessary jargon

## Site Architecture

### Learning Tracks
1. **Beginner Track**: "Getting Started"
   - Claude Code basics + simple bioinformatics tasks
   - File format handling and basic QC workflows
   - Understanding error messages and debugging

2. **Intermediate Track**: "Workflow Automation" 
   - Custom slash commands for common analyses
   - Project organization and CLAUDE.md setup
   - Pipeline integration and tool chaining

3. **Advanced Track**: "Power User Techniques"
   - Complex multi-omics workflows
   - Team collaboration and standardization
   - MCP integration and automation patterns

### Core Content Areas
- **Getting Started Guides**: Installation, setup, first analyses
- **Custom Commands Library**: Bioinformatics-specific slash commands
- **Workflow Examples**: RNA-seq, variant calling, single-cell, etc.
- **Project Templates**: CLAUDE.md examples and directory structures
- **Best Practices**: Code review, reproducibility, team standards
- **Troubleshooting**: Common issues and debugging strategies

## Technical Requirements

### Jekyll Configuration
- **Theme**: Documentation-focused theme with good code highlighting
- **Plugins**: Search, navigation, syntax highlighting, table of contents
- **Responsive Design**: Mobile-friendly for reference use
- **Performance**: Fast loading with optimized images and assets

### Custom Features
- **Copy-paste code blocks** with syntax highlighting
- **Download links** for templates and configuration files
- **Interactive command builder** for complex workflows
- **Progress tracking** for learning paths
- **Search functionality** across all content

## Team Collaboration

### Contributing Guidelines
- **Issue tracking** for content requests and bug reports
- **Pull request template** for content contributions
- **Review process** ensuring accuracy and educational value
- **Style guide** for consistent writing and formatting

### Community Engagement
- **Discussion forum** for questions and sharing workflows
- **Example submissions** from community members
- **Regular content updates** based on Claude Code feature releases
- **Feedback collection** for continuous improvement

## Development Phases

### Phase 1: Foundation (Current)
- [x] Project setup and planning documents
- [x] Git repository initialization and remote connection
- [ ] Basic Jekyll site structure
- [ ] Core planning documents completed

### Phase 2: Content Framework (Weeks 1-2)
- [ ] Educational progression designed
- [ ] Site navigation and structure implemented
- [ ] Template system for consistent content creation
- [ ] First beginner content modules

### Phase 3: Core Content Development (Weeks 3-8)
- [ ] Complete beginner track content
- [ ] Intermediate workflow examples
- [ ] Custom commands library
- [ ] Advanced integration patterns

### Phase 4: Community Features (Weeks 9-12)
- [ ] Interactive elements and tools
- [ ] Community contribution system
- [ ] Search and navigation optimization
- [ ] Mobile responsiveness and performance

### Phase 5: Launch & Iteration (Weeks 13+)
- [ ] Public launch and community outreach
- [ ] Feedback collection and content refinement
- [ ] Regular updates and new content additions
- [ ] Integration with Claude Code feature releases

## Success Metrics

### Educational Effectiveness
- **Clear learning progression** from novice to expert
- **Measurable skill development** through practical exercises
- **High user satisfaction** based on community feedback
- **Adoption by research teams** as standard reference

### Technical Quality
- **Fast site performance** (<3 second page loads)
- **High search visibility** for bioinformatics + Claude Code queries
- **Mobile compatibility** for reference use
- **Accessible design** meeting WCAG standards

### Community Impact
- **Active community engagement** through issues and discussions
- **Content contributions** from bioinformatics community
- **Citation and reference** in academic and industry contexts
- **Integration into training programs** at universities and institutes

## Important Files and References

### Planning Documents
- `DEVELOPMENT_ROADMAP.md` - Detailed project timeline and milestones
- `CONTENT_STRUCTURE.md` - Complete site architecture and navigation
- `EDUCATIONAL_FRAMEWORK.md` - Learning progression and pedagogical approach
- `TECHNICAL_REQUIREMENTS.md` - Jekyll setup and deployment specifications

### Source Content
- `bioinformatics-context-reference-guide.md` - Technical reference material
- `claude-code-best-practices.md` - Claude Code usage patterns
- `bioinformatics-one-liners.md` - Command-line examples and utilities

### Development Resources
- `.gitignore` - Git ignore patterns for Jekyll development
- `docs/` - GitHub Pages source content
- `assets/` - Images, CSS, and other static resources

## Current Status
**Phase 1: Foundation** - In Progress
- Git repository initialized and connected to remote
- Core planning documents in development
- Ready to begin Jekyll site structure and content framework

---

*This CLAUDE.md file serves as the central project documentation and should be updated as the project evolves. All team members should refer to this file for project context, workflow instructions, and current status.*