# Test_Setup/python/conftest.py

import pytest
from pathlib import Path
import pandas as pd

# Path to the CSV that the C++ sim writes
CSV_PATH = Path(__file__).parent.parent / "cpp" / "sim_results.csv"

@pytest.fixture(scope="module")
def df():
    """Load the simulation results into a pandas DataFrame."""
    assert CSV_PATH.exists(), f"CSV not found at: {CSV_PATH}"
    # read everything as floats where possible
    df = pd.read_csv(CSV_PATH)
    return df
