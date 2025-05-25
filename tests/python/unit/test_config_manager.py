"""
Unit tests for ConfigManager class.

Tests the configuration loading, validation, and management functionality
of the stellar plotter system.
"""

import pytest
import tempfile
import yaml
from pathlib import Path
from unittest.mock import patch, mock_open

# Import the module under test
from stellar_plotter import ConfigManager, PlotConfig


class TestConfigManager:
    """Test suite for ConfigManager class."""
    
    def test_init_with_no_config_file(self):
        """Test ConfigManager initialization with no config file."""
        config_manager = ConfigManager()
        assert config_manager.config is not None
        assert isinstance(config_manager.config, PlotConfig)
    
    def test_init_with_nonexistent_config_file(self):
        """Test ConfigManager initialization with non-existent config file."""
        config_manager = ConfigManager("nonexistent_file.yaml")
        assert config_manager.config is not None
        assert isinstance(config_manager.config, PlotConfig)
    
    def test_get_config_returns_plot_config(self):
        """Test that get_config returns a PlotConfig instance."""
        config_manager = ConfigManager()
        config = config_manager.get_config()
        
        assert isinstance(config, PlotConfig)
        assert hasattr(config, 'data_directory')
        assert hasattr(config, 'output_directory')
        assert hasattr(config, 'style')
        assert hasattr(config, 'validation')
        assert hasattr(config, 'performance')
    
    def test_validate_paths_with_existing_data_dir(self, temp_output_dir, sample_data_dir):
        """Test path validation with existing data directory."""
        # Create the sample data directory
        sample_data_dir.mkdir(parents=True, exist_ok=True)
        
        config_manager = ConfigManager()
        config = config_manager.get_config()
        config.data_directory = str(sample_data_dir)
        config.output_directory = str(temp_output_dir)
        
        # Update the config manager's config
        config_manager.config = config
        
        result = config_manager.validate_paths()
        assert result is True
        assert temp_output_dir.exists()
    
    def test_validate_paths_with_nonexistent_data_dir(self, temp_output_dir):
        """Test path validation with non-existent data directory."""
        config_manager = ConfigManager()
        config = config_manager.get_config()
        config.data_directory = "/nonexistent/directory"
        config.output_directory = str(temp_output_dir)
        
        # Update the config manager's config
        config_manager.config = config
        
        result = config_manager.validate_paths()
        assert result is False
    
    def test_validate_paths_creates_output_directory(self, sample_data_dir, temp_output_dir):
        """Test that validate_paths creates output directory if it doesn't exist."""
        # Create the sample data directory
        sample_data_dir.mkdir(parents=True, exist_ok=True)
        
        # Use a subdirectory that doesn't exist yet
        new_output_dir = temp_output_dir / "new_subdir"
        assert not new_output_dir.exists()
        
        config_manager = ConfigManager()
        config = config_manager.get_config()
        config.data_directory = str(sample_data_dir)
        config.output_directory = str(new_output_dir)
        
        # Update the config manager's config
        config_manager.config = config
        
        result = config_manager.validate_paths()
        assert result is True
        assert new_output_dir.exists()
    
    def test_config_manager_logger_setup(self):
        """Test that ConfigManager sets up logging correctly."""
        config_manager = ConfigManager()
        assert hasattr(config_manager, 'logger')
        assert config_manager.logger.name == 'stellar_plotter'


@pytest.mark.unit
class TestPlotConfigDefaults:
    """Test suite for PlotConfig default values."""
    
    def test_plot_config_defaults(self):
        """Test that PlotConfig has correct default values."""
        config = PlotConfig()
        
        # Test directory defaults
        assert config.data_directory == "../data/"
        assert config.output_directory == "../plots/"
        
        # Test style defaults
        assert config.style.theme == 'publication'
        assert config.style.figure_size == (10, 6)
        assert config.style.dpi == 300
        assert config.style.grid is True
        
        # Test validation defaults
        assert config.validation.min_mass_solar == 0.01
        assert config.validation.max_mass_solar == 10.0
        
        # Test performance defaults
        assert config.performance.enable_parallel_loading is True
        assert config.performance.chunk_size == 10
    
    def test_plot_config_eos_patterns(self):
        """Test that PlotConfig has correct EOS patterns."""
        config = PlotConfig()
        
        expected_patterns = {
            'hybrid': 'tov_solution_magnetic_bps_bbp_polytrope_*.csv',
            'neutron_relativistic': 'neutron_relativistic_rhoc_*.csv',
            'electron_relativistic': 'electron_relativistic_rhoc_*.csv',
            'electron_non_relativistic': 'electron_non_relativistic_rhoc_*.csv',
            'neutron_non_relativistic': 'neutron_non_relativistic_rhoc_*.csv'
        }
        
        assert config.eos_patterns == expected_patterns
    
    def test_cpp_encoded_density_parsing(self):
        """Test C++ encoded density parsing functionality."""
        config = PlotConfig()
        
        # Test positive exponent encoding
        assert config.parse_cpp_encoded_density("1.00pp16") == 1.00e+16
        assert config.parse_cpp_encoded_density("5.00pp18") == 5.00e+18
        
        # Test single 'p' for 'e' replacement
        assert config.parse_cpp_encoded_density("1.50p-03") == 1.50e-03
        
        # Test edge cases
        assert config.parse_cpp_encoded_density("1.0p15") == 1.0e15 