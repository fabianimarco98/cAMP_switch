# Contributing to cAMP_switch

Thank you for your interest in contributing to the cAMP binder design project!

## How to Contribute

### Reporting Issues

If you encounter bugs or have suggestions:
1. Check if the issue already exists
2. Create a new issue with:
   - Clear description of the problem
   - Steps to reproduce
   - Expected vs. actual behavior
   - System information (OS, software versions)
   - Error messages or logs

### Suggesting Enhancements

For feature requests:
1. Describe the proposed feature
2. Explain the use case and benefits
3. Provide examples if applicable

### Contributing Code

1. **Fork the repository**
   ```bash
   git clone https://github.com/your-username/cAMP_switch.git
   cd cAMP_switch
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make changes**
   - Follow existing code style
   - Add comments where necessary
   - Update documentation

4. **Test your changes**
   - Ensure scripts run without errors
   - Test on sample data if available
   - Verify documentation is accurate

5. **Commit changes**
   ```bash
   git add .
   git commit -m "Add feature: description"
   ```

6. **Push and create pull request**
   ```bash
   git push origin feature/your-feature-name
   ```

## Code Style Guidelines

### Python
- Follow PEP 8 style guide
- Use meaningful variable names
- Add docstrings to functions
- Include type hints where appropriate

Example:
```python
def analyze_structure(pdb_file: str, cutoff: float = 5.0) -> dict:
    """
    Analyze protein structure from PDB file.
    
    Args:
        pdb_file: Path to PDB file
        cutoff: Distance cutoff in Angstroms
    
    Returns:
        Dictionary with analysis results
    """
    # Implementation
    pass
```

### Shell Scripts
- Use bash shebang: `#!/bin/bash`
- Set error handling: `set -e`
- Add comments for complex commands
- Use meaningful variable names (uppercase for globals)

### Documentation
- Use Markdown for documentation
- Include examples and usage
- Keep documentation up-to-date with code changes
- Use clear, concise language

## Project Structure

When adding new features:
- Scripts go in appropriate `scripts/` directories
- Documentation in `docs/`
- Examples in `inputs/` directories
- Update README.md if adding major features

## Areas for Contribution

We welcome contributions in:

1. **Pipeline Improvements**
   - Optimization of existing scripts
   - Additional analysis tools
   - Automated workflows

2. **Documentation**
   - Tutorial improvements
   - Example workflows
   - Troubleshooting guides
   - Video tutorials

3. **Testing**
   - Unit tests for scripts
   - Integration tests
   - Test data sets

4. **Visualization**
   - Improved plotting scripts
   - Interactive visualizations
   - 3D structure viewers

5. **Performance**
   - Code optimization
   - Parallelization
   - Memory efficiency

## Development Setup

1. Install development dependencies:
```bash
pip install -r requirements.txt
pip install pytest black flake8  # For testing and linting
```

2. Run tests (when available):
```bash
pytest tests/
```

3. Check code style:
```bash
black --check .
flake8 .
```

## Pull Request Process

1. Update documentation for any user-facing changes
2. Add yourself to CONTRIBUTORS.md
3. Ensure all tests pass
4. Request review from maintainers
5. Address review comments
6. Squash commits if requested

## Questions?

Feel free to:
- Open an issue for discussion
- Contact maintainers
- Join community discussions

## Code of Conduct

Be respectful and constructive:
- Use welcoming and inclusive language
- Be respectful of differing viewpoints
- Accept constructive criticism gracefully
- Focus on what's best for the community

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Recognition

Contributors will be acknowledged in:
- CONTRIBUTORS.md file
- Release notes
- Project documentation

Thank you for contributing! 🎉
