"""
Unit tests for EOSDataProcessor class.

Tests the CSV file discovery, parsing, data validation, and processing
functionality of the stellar plotter system.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
import tempfile
import os
import sys
sys.path.append('../../scripts')

# Import the module under test
from stellar_plotter import EOSDataProcessor, EOSType, StellarModel, EOSDataset, PlotConfig


@pytest.mark.unit
class TestEOSDataProcessor:
    """Test suite for EOSDataProcessor class."""

    def test_init_with_config(self, mock_config):
        """Test EOSDataProcessor initialization with configuration."""
        # Convert mock_config dict to PlotConfig object
        config = PlotConfig()
        config.data_directory = mock_config.get('data_directory', '../data/')
        config.output_directory = mock_config.get('output_directory', '../plots/')

        processor = EOSDataProcessor(config)
        assert processor.config == config
        assert hasattr(processor, 'logger')
        assert processor.logger.name == 'stellar_plotter'

    def test_discover_files_hybrid_eos(self, multiple_hybrid_files, mock_config):
        """Test file discovery for hybrid EOS files."""
        # Convert mock_config dict to PlotConfig object and update data directory
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        processor = EOSDataProcessor(config)

        # discover_files expects list of strings, not EOSType enum
        files_dict = processor.discover_files(['hybrid'])

        assert 'hybrid' in files_dict
        assert len(files_dict['hybrid']) == 5  # We created 5 hybrid files
        assert all('tov_solution_magnetic_bps_bbp_polytrope_' in f for f in files_dict['hybrid'])
        assert all(f.endswith('.csv') for f in files_dict['hybrid'])

    def test_discover_files_neutron_relativistic(self, sample_neutron_csv_file, mock_config):
        """Test file discovery for neutron relativistic EOS files."""
        # Convert mock_config dict to PlotConfig object and update data directory
        config = PlotConfig()
        config.data_directory = str(sample_neutron_csv_file.parent)
        processor = EOSDataProcessor(config)

        files_dict = processor.discover_files(['neutron_relativistic'])

        assert 'neutron_relativistic' in files_dict
        assert len(files_dict['neutron_relativistic']) == 1
        assert 'neutron_relativistic_rhoc_5.00pp18.csv' in files_dict['neutron_relativistic'][0]

    def test_discover_files_empty_directory(self, temp_output_dir, mock_config):
        """Test file discovery in empty directory."""
        config = PlotConfig()
        config.data_directory = str(temp_output_dir)
        processor = EOSDataProcessor(config)

        files_dict = processor.discover_files(['hybrid'])

        # Should return empty dict or empty list for hybrid
        assert files_dict.get('hybrid', []) == []

    def test_discover_files_nonexistent_directory(self, mock_config):
        """Test file discovery in non-existent directory."""
        config = PlotConfig()
        config.data_directory = "/nonexistent/directory"
        processor = EOSDataProcessor(config)

        files_dict = processor.discover_files(['hybrid'])

        # Should return empty dict or empty list for hybrid
        assert files_dict.get('hybrid', []) == []

    def test_parse_central_density_hybrid(self, mock_config):
        """Test density parsing from hybrid EOS filename."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        # Test standard scientific notation
        density = processor.parse_central_density("tov_solution_magnetic_bps_bbp_polytrope_1.50e+16.csv", "hybrid")
        assert density == 1.50e+16

        # Test different density values
        density = processor.parse_central_density("tov_solution_magnetic_bps_bbp_polytrope_5.00e+15.csv", "hybrid")
        assert density == 5.00e+15

    def test_parse_central_density_neutron(self, mock_config):
        """Test density parsing from neutron relativistic filename with C++ encoding."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        # Test C++ encoded density (pp notation)
        density = processor.parse_central_density("neutron_relativistic_rhoc_1.00pp16.csv", "neutron_relativistic")
        assert density == 1.00e+16

        # Test different encoded values
        density = processor.parse_central_density("neutron_relativistic_rhoc_5.00pp18.csv", "neutron_relativistic")
        assert density == 5.00e+18

    def test_parse_central_density_invalid(self, mock_config):
        """Test density parsing from invalid filename."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        # Test filename without density
        density = processor.parse_central_density("invalid_filename.csv", "hybrid")
        assert density is None

        # Test filename with malformed density
        density = processor.parse_central_density("file_with_invalid_density_abc.csv", "hybrid")
        assert density is None

    def test_load_stellar_model_valid(self, sample_hybrid_csv_file, mock_config):
        """Test loading a valid stellar model."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        stellar_model = processor.load_stellar_model(str(sample_hybrid_csv_file), "hybrid", load_profile=False)

        assert isinstance(stellar_model, StellarModel)
        assert stellar_model.filename == sample_hybrid_csv_file.name
        assert stellar_model.eos_type == "hybrid"
        assert stellar_model.mass_solar > 0
        assert stellar_model.radius_km > 0
        assert stellar_model.central_density > 0
        assert stellar_model.log_central_density > 0

    def test_load_stellar_model_with_profile(self, sample_hybrid_csv_file, mock_config):
        """Test loading stellar model with profile data."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        stellar_model = processor.load_stellar_model(str(sample_hybrid_csv_file), "hybrid", load_profile=True)

        assert isinstance(stellar_model, StellarModel)
        assert stellar_model.log_r is not None
        assert stellar_model.log_m is not None
        assert stellar_model.log_p is not None
        assert len(stellar_model.log_r) == 100  # Sample data has 100 points

    def test_load_stellar_model_missing_file(self, mock_config):
        """Test loading from non-existent file."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        stellar_model = processor.load_stellar_model("nonexistent_file.csv", "hybrid")

        assert stellar_model is None

    def test_load_eos_dataset_single_type(self, multiple_hybrid_files, mock_config):
        """Test loading complete EOS dataset for single type."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        processor = EOSDataProcessor(config)

        dataset = processor.load_eos_dataset("hybrid")

        assert isinstance(dataset, EOSDataset)
        assert dataset.eos_type == "hybrid"
        assert len(dataset.models) == 5  # 5 hybrid files
        assert all(isinstance(model, StellarModel) for model in dataset.models)

    def test_load_eos_dataset_no_files(self, temp_output_dir, mock_config):
        """Test loading EOS dataset when no files exist."""
        config = PlotConfig()
        config.data_directory = str(temp_output_dir)
        processor = EOSDataProcessor(config)

        dataset = processor.load_eos_dataset("hybrid")

        assert dataset is None


@pytest.mark.unit
class TestEOSDataProcessorEdgeCases:
    """Test suite for EOSDataProcessor edge cases and error handling."""

    def test_empty_csv_file(self, sample_data_dir, mock_config):
        """Test handling of empty CSV file."""
        sample_data_dir.mkdir(parents=True, exist_ok=True)
        empty_file = sample_data_dir / "tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv"
        empty_file.write_text("")  # Empty file

        config = PlotConfig()
        processor = EOSDataProcessor(config)
        stellar_model = processor.load_stellar_model(str(empty_file), "hybrid")

        assert stellar_model is None

        # Cleanup
        empty_file.unlink()

    def test_csv_file_with_nan_values(self, sample_data_dir, mock_config):
        """Test handling of CSV file with NaN values."""
        sample_data_dir.mkdir(parents=True, exist_ok=True)
        nan_file = sample_data_dir / "tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv"

        # Create data with NaN values
        data_with_nan = pd.DataFrame({
            'log_r[cm]': [5.0, np.nan, 6.0],
            'log_m[g]': [33.0, 33.2, np.nan],
            'log_P[dyne/cm^2]': [35.0, 30.0, 25.0]
        })
        data_with_nan.to_csv(nan_file, index=False)

        config = PlotConfig()
        processor = EOSDataProcessor(config)
        stellar_model = processor.load_stellar_model(str(nan_file), "hybrid")

        # Should handle NaN values gracefully - either load successfully or return None
        if stellar_model is not None:
            # If loaded, should have valid data
            assert stellar_model.mass_solar > 0
            assert stellar_model.radius_km > 0

        # Cleanup
        nan_file.unlink()

    def test_csv_file_missing_columns(self, sample_data_dir, mock_config):
        """Test handling of CSV file with missing required columns."""
        sample_data_dir.mkdir(parents=True, exist_ok=True)
        invalid_file = sample_data_dir / "tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv"

        # Create data with wrong column names
        invalid_data = pd.DataFrame({
            'wrong_column1': [1, 2, 3],
            'wrong_column2': [4, 5, 6]
        })
        invalid_data.to_csv(invalid_file, index=False)

        config = PlotConfig()
        processor = EOSDataProcessor(config)
        stellar_model = processor.load_stellar_model(str(invalid_file), "hybrid")

        # Should return None for invalid files
        assert stellar_model is None

        # Cleanup
        invalid_file.unlink()

    def test_very_large_dataset(self, mock_config):
        """Test handling of very large dataset (memory management)."""
        config = PlotConfig()
        processor = EOSDataProcessor(config)

        # Create a large synthetic stellar model with correct constructor
        stellar_model = StellarModel(
            mass_solar=2.0,
            radius_km=12.0,
            central_density=1e16,
            log_central_density=16.0,
            eos_type="hybrid",
            filename="large_test.csv"
        )

        # Add large profile data
        stellar_model.log_r = np.linspace(3, 6, 10000)
        stellar_model.log_m = np.linspace(30, 34, 10000)
        stellar_model.log_p = np.linspace(35, 20, 10000)

        # Should handle large datasets without issues - test validation
        if config.validation:
            is_valid, warnings = config.validation.validate_model(stellar_model)
            assert isinstance(is_valid, bool)  # Should complete without memory errors
