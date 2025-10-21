# Repository Guidelines

## Project Structure & Module Organization
The core Python package lives in `nisarhdf/`, with product-specific readers (`nisarRSLCHDF.py`, `nisarGUNWHDF.py`, etc.) sharing the base abstractions in `nisarBaseHDF.py`. Command-line entry points (`nisarH5ToImage.py`, `multilookRSLC.py`, `s3Listing.py`) expose the functionality described in `pyproject.toml`. Support notebooks for discovery and demos are kept in `Notebooks/`, and example automation lives in `scripts/nisar_gcov_cog.py`. Environment setup is captured in `environment.yml`; distribution metadata sits in `pyproject.toml`.

## Build, Test, and Development Commands
Create a virtual env before editing: `python -m venv .venv && source .venv/bin/activate`. Install editable deps with `pip install -e .`. For ad-hoc builds, `python -m build` emits a wheel. Run the main converter via `nisarh5toimage --help`, or `python -m nisarhdf.nisarH5ToImage` when debugging modules. `scripts/nisar_gcov_cog.py` expects GDAL bindings; validate availability with `gdalinfo --version` inside your env.

## Coding Style & Naming Conventions
Match the existing modules: 4-space indentation, snake_case functions, and CapWords classes (`NisarRSLCHDF`). Keep imports grouped stdlib/third-party/local, and prefer explicit exports via `__all__` when adding modules. Write docstrings describing I/O and data shapes, and surface logging through `logging` instead of new print statements. Preserve NumPy vectorization and avoid per-pixel Python loops unless profiling shows gains.

## Testing Guidelines
There is no automated test suite yet; prefer `pytest` and name new files `test_<module>.py` under a `tests/` package. When touching readers, add minimal regression fixtures using representative `.h5` files and assert on key arrays (e.g., `unwrappedPhase.shape`). For CLI work, exercise `nisarh5toimage` against a lightweight sample and capture the JSON summary for review. Document manual validation steps in pull requests until coverage thresholds are defined.

## Commit & Pull Request Guidelines
Recent commits are single-line, present-tense summaries (e.g., `add multilookRSLC handling`). Follow that style and keep body text for context, including data sources touched. Every pull request should state the problem, summarize the solution, and reference related issues or tickets. Attach before/after metrics or imagery when altering performance-sensitive code, and note how you validated the change (command snippets, datasets). Request review from a domain maintainer before merging.

## Data & Environment Notes
Large NISAR sample files are not tracked; store paths in `.gitignore` and document download commands instead. Use `environment.yml` when replicating the full geospatial stack (`conda env create -f environment.yml`) to ensure GDAL/rasterio compatibility. Confirm AWS credentials and other secrets are redacted before sharing notebooks or logs.
