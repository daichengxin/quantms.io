# Contributing to QPX

Thank you for your interest in contributing to QPX!

## Code of Conduct

We follow the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

## Getting Started

### Setup

1. Fork and clone the repository:

   ```bash
   git clone https://github.com/your-username/qpx.git
   cd qpx
   git remote add upstream https://github.com/bigbio/qpx.git
   ```

2. Install dependencies (Python 3.10+):

   ```bash
   pip install -e ".[dev]"
   ```

### Development Workflow

Create a feature branch:

```bash
git checkout -b feature/your-feature-name
```

## Code Standards

Format code with ruff:

```bash
ruff check .
ruff format .
```

## Testing

Run tests (excludes integration tests by default):

```bash
pytest -vv -m "not integration"
```

Test markers:

- `integration`: Integration tests requiring external resources
- `large_data`: Tests requiring large data files

## Submitting Changes

1. Commit and push your changes:

   ```bash
   git add .
   git commit -m "Brief description of changes"
   git push origin feature/your-feature-name
   ```

2. Create a pull request on GitHub with:
   - Clear description of changes
   - Related issue numbers (e.g., "Fixes #123")
   - Test results

## Reporting Issues

**Bugs**: Include title, description, steps to reproduce, environment (QPX version, Python version, OS), and error messages.

**Features**: Describe the use case, proposed solution, and alternatives considered.

## Documentation

Build and preview documentation locally:

```bash
mkdocs serve  # Available at http://127.0.0.1:8000/
```

## Contact

- **Maintainers**: Yasset Perez-Riverol (PRIDE Team, EMBL-EBI), Ping Zheng (Chongqing University)
- **Issues**: [GitHub Issues](https://github.com/bigbio/qpx/issues)

## License

Contributions are licensed under [Apache License 2.0](LICENSE). Please add your name to the contributors list in README.md.
