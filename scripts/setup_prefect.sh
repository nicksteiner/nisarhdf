#!/bin/bash
# Quick start: install Prefect and run the example flow

set -e

echo "=== Installing Prefect ==="
pip install -q prefect

echo ""
echo "=== Testing standalone logging ==="
python scripts/test_logging.py 2>&1 | grep nisar_gcov_cog

echo ""
echo "=== Prefect Flow Example ==="
echo "Edit scripts/prefect_gcov_flow.py to set your input directory, then run:"
echo ""
echo "  python scripts/prefect_gcov_flow.py"
echo ""
echo "Or start Prefect server for UI:"
echo ""
echo "  prefect server start       # http://localhost:4200"
echo "  python scripts/prefect_gcov_flow.py"
echo ""
echo "All logs from nisar_gcov_cog.py will appear in Prefect UI!"
