# Contributing

- Write READMEs, docstrings and comments in English.
- Follow PEP 8 and give public modules, classes and functions a docstring.
- Raise and catch specific exceptions; do not use a bare `except`.
- Report diagnostics with the `logging` module rather than `print`.
- Add pytest tests for new features and bug fixes.
- Use comments to explain non-obvious logic, scientific assumptions and units.
  Describe the current behaviour, not how it changed; put that in the commit
  message.
- Declare Python dependencies in the package's `pyproject.toml`. Document
  external programs (Gaussian, AmberTools, GROMACS) in the README.

Run the tests from the respective directory:

```bash
cd electrofit-workspace && python -m pytest                          # electrofit
cd electrofit-workspace/packages/ip6-analysis && python -m pytest    # ip6-analysis
```
