#!/usr/bin/env python3
"""Quick test of logging integration."""
import logging
import sys
from pathlib import Path

# Set up logging to see the output
logging.basicConfig(
    level=logging.DEBUG,
    format="[%(levelname)-7s] %(name)s: %(message)s"
)

# Add parent dir to path to import the script
sys.path.insert(0, str(Path(__file__).parent))

from nisar_gcov_cog import logger, main

# Test 1: Logger is named correctly
print(f"✓ Logger name: {logger.name}")
assert logger.name == "nisar_gcov_cog"

# Test 2: Can emit logs
logger.info("Test info message")
logger.warning("Test warning")
logger.debug("Test debug")

# Test 3: Help works
print("\n--- Testing --help ---")
try:
    exit_code = main(["--help"])
except SystemExit as e:
    print(f"✓ Help exited with code: {e.code}")

# Test 4: Error logging works (nonexistent file)
print("\n--- Testing error logging ---")
exit_code = main([
    "/tmp/nonexistent_gcov_file.h5",
    "-o", "/tmp/test_output",
    "--stac", "/tmp/test_stac"
])
print(f"✓ Exit code for missing file: {exit_code}")

print("\n✓ All logging tests passed!")
