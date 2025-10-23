# Summary: Prefect Integration for nisar_gcov_cog.py

## Changes Made

### 1. **Replaced print() with logging**
   - Added `import logging` and module-level `logger = logging.getLogger(__name__)`
   - Replaced all `print()` → `logger.info()`, `logger.warning()`, `logger.error()`, `logger.debug()`
   - Uses lazy formatting: `logger.info("Processing %s", filename)` instead of f-strings

### 2. **Smart logging setup**
   - Added `setup_logging(verbose)` function that only configures if no handlers exist
   - When Prefect (or any framework) has already configured logging, the script respects that
   - Standalone CLI still gets nice console output automatically

### 3. **Updated __main__ block**
   ```python
   if __name__ == "__main__":
       args = parse_args(sys.argv[1:])
       setup_logging(verbose=args.verbose)  # Only when standalone
       raise SystemExit(main(sys.argv[1:]))
   ```

---

## How to Use with Prefect

### Option 1: Simple task wrapper (see `scripts/prefect_gcov_flow.py`)
```python
from prefect import task, get_run_logger
from scripts.nisar_gcov_cog import main

@task
def convert_gcov(h5_path, output_dir):
    logger = get_run_logger()
    logger.info("Converting %s", h5_path)
    exit_code = main([str(h5_path), "-o", str(output_dir), ...])
    return exit_code == 0
```

**Benefits:**
- All logs appear in Prefect UI with task context
- Automatic retry/failure tracking
- No double logging (Prefect's handler intercepts everything)

### Option 2: Direct import and refactor
If you want more granular tasks, import the helper functions:
```python
from scripts.nisar_gcov_cog import (
    open_frequency_readers,
    export_frequency_assets,
    build_item,
)

@task
def read_gcov(h5_path):
    return open_frequency_readers(h5_path, sigma0=False, db=True)

@task
def export_cogs(reader, freq, h5_path, output):
    return export_frequency_assets(reader, freq, h5_path, output)
```

---

## Testing

Run the test script:
```bash
cd scripts
python test_logging.py
```

Expected output shows structured logs with proper levels:
```
[INFO   ] nisar_gcov_cog: Test info message
[WARNING] nisar_gcov_cog: Test warning
[ERROR  ] nisar_gcov_cog: Input not found: /tmp/nonexistent.h5
```

---

## Files Created

1. **`scripts/nisar_gcov_cog.py`** — refactored with logging (edited in place)
2. **`scripts/prefect_gcov_flow.py`** — complete Prefect flow example with batch processing
3. **`scripts/LOGGING_INTEGRATION.md`** — comprehensive guide for all frameworks (Airflow, Dagster, etc.)
4. **`scripts/test_logging.py`** — quick smoke test for the logging setup

---

## Next Steps

1. **Try Prefect locally:**
   ```bash
   pip install prefect
   python scripts/prefect_gcov_flow.py
   ```

2. **Deploy to Prefect Cloud/Server:**
   - Logs automatically stream to Prefect UI
   - Set up schedules, notifications, and monitoring

3. **Adapt for your framework:**
   - See `LOGGING_INTEGRATION.md` for Airflow/Dagster patterns
   - Core principle: import `main()`, call with argv list, logs propagate

---

## Key Advantages

✅ **Zero breaking changes** — CLI still works exactly the same  
✅ **Framework-agnostic** — works with any Python logging consumer  
✅ **Structured** — proper log levels instead of print noise  
✅ **Traceable** — Prefect UI shows which task emitted which log  
✅ **Retryable** — Prefect can retry tasks with full log context  
✅ **Observable** — metrics, timing, and logs in one dashboard
