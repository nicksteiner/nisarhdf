# Logging Integration for nisar_gcov_cog.py

## Overview
`nisar_gcov_cog.py` now uses Python's standard `logging` module instead of `print()` statements, making it compatible with orchestration frameworks like **Prefect**, **Airflow**, **Dagster**, or any custom logging infrastructure.

---

## Key Changes

### 1. Module-level logger
```python
import logging
logger = logging.getLogger(__name__)
```
- Uses `__name__` so logs are namespaced as `nisar_gcov_cog`
- Propagates to parent loggers (e.g., Prefect's root logger)
- No side effects when imported as a library

### 2. Log levels
- `logger.info()` → progress messages (e.g., "Processing file.h5")
- `logger.warning()` → non-fatal issues (e.g., fallback to in-memory write)
- `logger.error()` → file not found, conversion failures
- `logger.debug()` → verbose output (conversion min/max stats)

### 3. Standalone CLI setup
When run directly (`python nisar_gcov_cog.py ...`), the script configures a basic console logger:
```python
def setup_logging(verbose: bool = False):
    if logging.getLogger().hasHandlers():
        return  # Already configured by a parent (e.g., Prefect)
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, ...)
```
If the root logger already has handlers (set by Prefect/Airflow/etc.), it skips setup and respects the parent configuration.

---

## Usage Examples

### Direct CLI (standalone)
```bash
# Standard logging (INFO level)
python ./scripts/nisar_gcov_cog.py input.h5 -o ./outputs --stac ./stac

# Verbose logging (DEBUG level)
python ./scripts/nisar_gcov_cog.py input.h5 -o ./outputs --stac ./stac --verbose
```

### Prefect Integration
See `scripts/prefect_gcov_flow.py` for a complete example.

```python
from pathlib import Path
from prefect import flow, task, get_run_logger
from scripts.nisar_gcov_cog import main  # imports without side effects

@task
def convert_gcov(h5_path: Path, output_dir: Path):
    logger = get_run_logger()
    logger.info("Converting %s", h5_path.name)
    
    exit_code = main([str(h5_path), "-o", str(output_dir), "--stac", "./stac"])
    return exit_code == 0

@flow
def batch_convert(files):
    for f in files:
        convert_gcov(f, Path("./outputs"))
```

**Benefits:**
- Prefect automatically captures all `logger.*` calls
- Logs appear in Prefect UI with task context (retries, timing, etc.)
- No double-logging (Prefect's handler takes precedence)
- Structured logs with timestamps, levels, and context

---

## Running the Prefect Example

1. **Install Prefect** (if not already):
   ```bash
   pip install prefect
   ```

2. **Start Prefect server** (optional, for UI):
   ```bash
   prefect server start
   ```

3. **Run the flow**:
   ```bash
   python ./scripts/prefect_gcov_flow.py
   ```

4. **View logs**:
   - In terminal: structured logs with Prefect context
   - In Prefect UI (http://localhost:4200): task-level logs with filtering

---

## Airflow / Dagster / Custom Frameworks

The pattern is identical:
1. Import `nisar_gcov_cog.main` (don't invoke `setup_logging()`)
2. Configure your framework's root logger
3. Call `main(argv_list)` from a task/operator
4. All logs propagate to your framework's handlers

**Example (Airflow PythonOperator):**
```python
from airflow.operators.python import PythonOperator

def airflow_task(**context):
    from scripts.nisar_gcov_cog import main
    ti = context['ti']
    ti.log.info("Starting GCOV conversion")
    
    exit_code = main(["input.h5", "-o", "./outputs", "--stac", "./stac"])
    ti.log.info("Exit code: %d", exit_code)

task = PythonOperator(task_id="convert_gcov", python_callable=airflow_task)
```

---

## Configuration Tips

### Adjusting log format
When calling `setup_logging()` standalone, customize format:
```python
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s | %(name)s | %(message)s",  # simpler format
    handlers=[logging.StreamHandler(), logging.FileHandler("gcov.log")]  # file + console
)
```

### Silencing GDAL warnings
GDAL can be noisy; add before imports:
```python
import logging
logging.getLogger("osgeo").setLevel(logging.ERROR)
```

### JSON logs for cloud ingestion
Use `python-json-logger`:
```bash
pip install python-json-logger
```
```python
from pythonjsonlogger import jsonlogger

handler = logging.StreamHandler()
handler.setFormatter(jsonlogger.JsonFormatter())
logging.getLogger().addHandler(handler)
```

---

## Testing Log Integration

Quickly verify logs are structured:
```bash
python -c "
import logging
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(name)s: %(message)s')
from scripts.nisar_gcov_cog import main
main(['nonexistent.h5', '-o', '/tmp/test', '--stac', '/tmp/stac'])
"
```
Expected output:
```
[ERROR] nisar_gcov_cog: Input not found: nonexistent.h5
```

---

## Summary

✅ **Prefect-ready**: logs propagate without modification  
✅ **Standalone-friendly**: auto-configures when run directly  
✅ **Structured**: uses proper log levels (info/warn/error/debug)  
✅ **Framework-agnostic**: works with any Python logging consumer  
✅ **No breaking changes**: CLI behavior unchanged; just better logs
