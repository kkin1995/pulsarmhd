"""
Pytest configuration and fixtures for stellar plotter testing.

This module provides common fixtures and utilities for testing the
stellar structure plotting system.
"""

import pytest
import tempfile
import shutil
import os
import sys
from pathlib import Path
from typing import Dict, Any
import pandas as pd
import numpy as np

# Add scripts directory to Python path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

@pytest.fixture(scope="session")
def project_root():
    """Get the project root directory."""
    return Path(__file__).parent.parent.parent

@pytest.fixture(scope="session") 
def scripts_dir(project_root):
    """Get the scripts directory."""
    return project_root / "scripts"

@pytest.fixture
def temp_output_dir():
    """Provide temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp(prefix="stellar_plotter_test_")
    yield Path(temp_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)

@pytest.fixture
def sample_data_dir():
    """Provide sample CSV data directory for testing."""
    return Path(__file__).parent / "fixtures" / "sample_data"

@pytest.fixture
def mock_configs_dir():
    """Provide mock configuration files directory."""
    return Path(__file__).parent / "fixtures" / "mock_configs"

@pytest.fixture
def expected_outputs_dir():
    """Provide expected outputs directory for comparison."""
    return Path(__file__).parent / "fixtures" / "expected_outputs"

@pytest.fixture
def mock_config():
    """Provide mock configuration dictionary for testing."""
    return {
        'data_directory': str(Path(__file__).parent / "fixtures" / "sample_data"),
        'output_directory': 'temp_test_output',
        'style': {
            'theme': 'publication',
            'figure_size': [10, 6],
            'dpi': 300,
            'grid': True,
            'alpha': 0.7,
            'marker_size': 30.0,
            'line_width': 2.0
        },
        'validation': {
            'min_mass_solar': 0.01,
            'max_mass_solar': 10.0,
            'min_radius_km': 1.0,
            'max_radius_km': 100.0
        },
        'performance': {
            'enable_parallel_loading': True,
            'chunk_size': 10
        }
    }

@pytest.fixture
def sample_stellar_model_data():
    """Generate sample stellar model data for testing."""
    # Create realistic stellar structure data
    n_points = 100
    log_r = np.linspace(3, 6, n_points)  # log10(radius) from 1km to 1000km
    log_m = np.linspace(30, 33.3, n_points)  # log10(mass) realistic range
    log_p = np.linspace(35, 25, n_points)  # log10(pressure) decreasing outward
    
    return pd.DataFrame({
        'log_r[cm]': log_r,
        'log_m[g]': log_m,
        'log_P[dyne/cm^2]': log_p
    })

@pytest.fixture
def sample_hybrid_csv_file(sample_data_dir, sample_stellar_model_data):
    """Create a sample hybrid EOS CSV file for testing."""
    sample_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a realistic filename following C++ convention
    filename = "tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv"
    filepath = sample_data_dir / filename
    
    # Save the sample data
    sample_stellar_model_data.to_csv(filepath, index=False)
    
    yield filepath
    
    # Cleanup with error handling
    try:
        if filepath.exists():
            filepath.unlink()
    except (PermissionError, OSError):
        pass

@pytest.fixture
def sample_neutron_csv_file(sample_data_dir, sample_stellar_model_data):
    """Create a sample neutron relativistic EOS CSV file for testing."""
    sample_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create filename with C++ encoded density notation
    filename = "neutron_relativistic_rhoc_5.00pp18.csv"
    filepath = sample_data_dir / filename
    
    # Save the sample data
    sample_stellar_model_data.to_csv(filepath, index=False)
    
    yield filepath
    
    # Cleanup with error handling
    try:
        if filepath.exists():
            filepath.unlink()
    except (PermissionError, OSError):
        pass

@pytest.fixture
def multiple_hybrid_files(sample_data_dir, sample_stellar_model_data):
    """Create multiple hybrid EOS files for dataset testing."""
    sample_data_dir.mkdir(parents=True, exist_ok=True)
    
    densities = ["1.00e+15", "5.00e+15", "1.00e+16", "5.00e+16", "1.00e+17"]
    filepaths = []
    
    for density in densities:
        filename = f"tov_solution_magnetic_bps_bbp_polytrope_{density}.csv"
        filepath = sample_data_dir / filename
        
        # Vary the data slightly for each density
        data = sample_stellar_model_data.copy()
        # Add some realistic variation based on density
        density_factor = float(density.replace('e+', 'e'))
        data['log_m[g]'] *= (1 + np.log10(density_factor/1e15) * 0.1)
        data['log_r[cm]'] *= (1 - np.log10(density_factor/1e15) * 0.05)
        
        data.to_csv(filepath, index=False)
        filepaths.append(filepath)
    
    yield filepaths
    
    # Cleanup with error handling for WSL permission issues
    for filepath in filepaths:
        try:
            if filepath.exists():
                filepath.unlink()
        except (PermissionError, OSError):
            # Ignore permission errors in WSL environment
            pass

@pytest.fixture
def mock_yaml_config(mock_configs_dir, temp_output_dir):
    """Create a mock YAML configuration file for testing."""
    mock_configs_dir.mkdir(parents=True, exist_ok=True)
    
    config_content = f"""
# Test configuration for stellar plotter
data_directory: "{Path(__file__).parent / 'fixtures' / 'sample_data'}"
output_directory: "{temp_output_dir}"

style:
  theme: "publication"
  figure_size: [10, 6]
  dpi: 300
  grid: true
  alpha: 0.7
  marker_size: 30.0
  line_width: 2.0

validation:
  min_mass_solar: 0.01
  max_mass_solar: 10.0
  min_radius_km: 1.0
  max_radius_km: 100.0

performance:
  enable_parallel_loading: true
  chunk_size: 10
"""
    
    config_file = mock_configs_dir / "test_config.yaml"
    config_file.write_text(config_content)
    
    yield config_file
    
    # Cleanup with error handling
    try:
        if config_file.exists():
            config_file.unlink()
    except (PermissionError, OSError):
        pass

@pytest.fixture(autouse=True)
def setup_test_environment(monkeypatch, temp_output_dir):
    """Set up test environment with proper paths and cleanup."""
    # Ensure matplotlib doesn't try to display plots during testing
    monkeypatch.setenv("MPLBACKEND", "Agg")
    
    # Set up temporary directories
    os.makedirs(temp_output_dir, exist_ok=True)
    
    yield
    
    # Cleanup is handled by temp_output_dir fixture

# Custom markers for test organization
def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line("markers", "unit: Unit tests for individual components")
    config.addinivalue_line("markers", "integration: Integration tests between components")
    config.addinivalue_line("markers", "performance: Performance and benchmark tests")
    config.addinivalue_line("markers", "visual: Visual regression tests for plots")
    config.addinivalue_line("markers", "slow: Tests that take longer to run")

def pytest_collection_modifyitems(config, items):
    """Automatically mark slow tests."""
    for item in items:
        # Mark tests with 'benchmark' in name as slow
        if "benchmark" in item.name or "performance" in item.name:
            item.add_marker(pytest.mark.slow)
        
        # Mark visual tests
        if "visual" in item.name or "plot" in item.name:
            item.add_marker(pytest.mark.visual) 