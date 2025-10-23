import os
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--gcov", action="store", default=None,
        help="Path to a real NISAR GCOV .h5 file (or set NISAR_GCOV_PATH)"
    )


@pytest.fixture(scope="session")
def gcov_path(request):
    # Priority: CLI option > env var
    cli = request.config.getoption("--gcov")
    env = os.getenv("NISAR_GCOV_PATH")
    path = cli or env
    if not path:
        pytest.skip("No GCOV path provided. Set NISAR_GCOV_PATH or use --gcov.")
    if not (path.startswith("s3://") or path.startswith("https://") or os.path.exists(path)):
        pytest.skip(f"GCOV path not accessible: {path}")
    return path
