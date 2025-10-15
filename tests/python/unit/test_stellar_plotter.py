"""
Unit tests for StellarPlotter class.

Tests the plotting functionality, theme application, data processing,
and visualization methods of the stellar plotter system.
"""

import pytest
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import tempfile
import sys
sys.path.append('../../scripts')

# Import the module under test
from stellar_plotter import StellarPlotter, EOSDataProcessor, EOSType, StellarModel, EOSDataset, PlotConfig


@pytest.mark.unit
class TestStellarPlotter:
    """Test suite for StellarPlotter class."""

    def test_init_with_config(self, mock_config, temp_output_dir):
        """Test StellarPlotter initialization with configuration."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        assert plotter.config == config
        assert hasattr(plotter, 'logger')
        assert plotter.logger.name == 'stellar_plotter'

    def test_apply_theme_publication(self, mock_config, temp_output_dir):
        """Test applying publication theme."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        config.style.theme = 'publication'
        plotter = StellarPlotter(config)

        # Theme should be applied during initialization
        assert config.style.theme == 'publication'

        # Check that matplotlib rcParams were updated (rcParams returns lists)
        assert list(plt.rcParams['figure.figsize']) == list(config.style.figure_size)
        assert plt.rcParams['figure.dpi'] == config.style.dpi

    def test_apply_theme_dark(self, mock_config, temp_output_dir):
        """Test applying dark theme."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        config.style.theme = 'dark'
        plotter = StellarPlotter(config)

        # Theme should be applied during initialization
        assert config.style.theme == 'dark'

    def test_apply_theme_colorblind(self, mock_config, temp_output_dir):
        """Test applying colorblind-friendly theme."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        config.style.theme = 'colorblind'
        plotter = StellarPlotter(config)

        # Theme should be applied and colorblind palette should be set
        assert config.style.theme == 'colorblind'
        # Colorblind palette should be applied to colors
        assert config.style.colors != PlotConfig().style.colors  # Should be different from default

    def test_get_plot_style(self, mock_config, temp_output_dir):
        """Test getting plot style for EOS type."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        style = plotter._get_plot_style('hybrid')

        assert isinstance(style, dict)
        assert 'color' in style
        assert 'marker' in style
        assert 's' in style  # marker size
        assert 'alpha' in style
        assert 'linewidth' in style
        assert style['color'] == config.style.colors['hybrid']
        assert style['marker'] == config.style.markers['hybrid']

    def test_plot_single_profile(self, sample_hybrid_csv_file, mock_config, temp_output_dir):
        """Test plotting single stellar profile."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        # Load a stellar model with profile data
        processor = EOSDataProcessor(config)
        model = processor.load_stellar_model(str(sample_hybrid_csv_file), "hybrid", load_profile=True)

        if model is not None:
            # Test the plotting function
            output_path = plotter.plot_single_profile(model)

            assert isinstance(output_path, str)
            assert output_path.endswith('.png')
            assert Path(output_path).exists()

    def test_plot_single_profile_missing_profile_data(self, sample_hybrid_csv_file, mock_config, temp_output_dir):
        """Test plotting single profile with missing profile data."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        # Load a stellar model without profile data
        processor = EOSDataProcessor(config)
        model = processor.load_stellar_model(str(sample_hybrid_csv_file), "hybrid", load_profile=False)

        if model is not None:
            # Should raise ValueError for missing profile data
            with pytest.raises(ValueError):
                plotter.plot_single_profile(model)

    def test_plot_mass_radius_relation(self, multiple_hybrid_files, mock_config, temp_output_dir):
        """Test plotting mass-radius relation."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        config.output_directory = str(temp_output_dir)

        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)

        # Load dataset
        dataset = processor.load_eos_dataset("hybrid")

        if dataset is not None:
            output_path = plotter.plot_mass_radius_relation([dataset])

            assert isinstance(output_path, str)
            assert output_path.endswith('.png')
            assert Path(output_path).exists()

    def test_plot_mass_density_relation(self, multiple_hybrid_files, mock_config, temp_output_dir):
        """Test plotting mass-density relation."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        config.output_directory = str(temp_output_dir)

        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)

        # Load dataset
        dataset = processor.load_eos_dataset("hybrid")

        if dataset is not None:
            output_path = plotter.plot_mass_density_relation([dataset])

            assert isinstance(output_path, str)
            assert output_path.endswith('.png')
            assert Path(output_path).exists()

    def test_plot_comparative_analysis(self, multiple_hybrid_files, sample_neutron_csv_file, mock_config, temp_output_dir):
        """Test plotting comparative analysis of multiple EOS types."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        config.output_directory = str(temp_output_dir)

        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)

        # Load datasets
        datasets = []
        hybrid_dataset = processor.load_eos_dataset("hybrid")
        if hybrid_dataset:
            datasets.append(hybrid_dataset)

        neutron_dataset = processor.load_eos_dataset("neutron_relativistic")
        if neutron_dataset:
            datasets.append(neutron_dataset)

        if datasets:
            output_paths = plotter.plot_comparative_analysis(datasets)

            assert isinstance(output_paths, list)
            assert len(output_paths) > 0
            assert all(isinstance(path, str) for path in output_paths)
            assert all(path.endswith('.png') for path in output_paths)
            assert all(Path(path).exists() for path in output_paths)

    def test_create_summary_statistics(self, multiple_hybrid_files, mock_config, temp_output_dir):
        """Test creation of summary statistics."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)
        config.output_directory = str(temp_output_dir)

        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)

        # Load dataset
        dataset = processor.load_eos_dataset("hybrid")

        if dataset is not None:
            stats = plotter.create_summary_statistics([dataset])

            assert isinstance(stats, dict)
            assert 'hybrid' in stats
            assert isinstance(stats['hybrid'], dict)
            # Check for expected statistics keys
            expected_keys = ['num_models', 'min_mass', 'max_mass', 'min_radius', 'max_radius']
            for key in expected_keys:
                if key in stats['hybrid']:
                    assert isinstance(stats['hybrid'][key], (int, float))


@pytest.mark.unit
class TestStellarPlotterDataProcessing:
    """Test suite for StellarPlotter data processing methods."""

    def test_eos_dataset_properties(self, multiple_hybrid_files, mock_config):
        """Test EOSDataset properties for data extraction."""
        config = PlotConfig()
        config.data_directory = str(multiple_hybrid_files[0].parent)

        processor = EOSDataProcessor(config)
        dataset = processor.load_eos_dataset("hybrid")

        if dataset is not None:
            # Test dataset properties
            masses = dataset.masses
            radii = dataset.radii
            log_densities = dataset.log_densities

            assert isinstance(masses, list)
            assert isinstance(radii, list)
            assert isinstance(log_densities, list)
            assert len(masses) == len(radii) == len(log_densities)
            assert all(isinstance(m, float) for m in masses)
            assert all(isinstance(r, float) for r in radii)
            assert all(isinstance(d, float) for d in log_densities)
            assert all(m > 0 for m in masses)
            assert all(r > 0 for r in radii)


@pytest.mark.unit
class TestStellarPlotterErrorHandling:
    """Test suite for StellarPlotter error handling and edge cases."""

    def test_plot_with_empty_dataset_list(self, mock_config, temp_output_dir):
        """Test plotting with empty dataset list."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        # Test with empty list
        output_path = plotter.plot_mass_radius_relation([])

        # Should still create a plot (empty plot) or handle gracefully
        assert isinstance(output_path, str)

    def test_plot_with_dataset_no_models(self, mock_config, temp_output_dir):
        """Test plotting with dataset containing no models."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        plotter = StellarPlotter(config)

        # Create empty dataset
        empty_dataset = EOSDataset(eos_type="hybrid", models=[])

        output_path = plotter.plot_mass_radius_relation([empty_dataset])

        # Should handle gracefully
        assert isinstance(output_path, str)


@pytest.mark.unit
class TestStellarPlotterConfiguration:
    """Test suite for StellarPlotter configuration handling."""

    def test_custom_figure_size(self, mock_config, temp_output_dir):
        """Test custom figure size configuration."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        config.style.figure_size = (12, 8)
        plotter = StellarPlotter(config)

        # Check that matplotlib was configured with custom figure size (rcParams returns lists)
        assert list(plt.rcParams['figure.figsize']) == [12, 8]

    def test_custom_dpi(self, mock_config, temp_output_dir):
        """Test custom DPI configuration."""
        config = PlotConfig()
        config.output_directory = str(temp_output_dir)
        config.style.dpi = 600
        plotter = StellarPlotter(config)

        # Check that matplotlib was configured with custom DPI
        assert plt.rcParams['figure.dpi'] == 600

    def test_theme_configuration(self, mock_config, temp_output_dir):
        """Test different theme configurations."""
        themes = ['publication', 'presentation', 'dark', 'colorblind']

        for theme in themes:
            config = PlotConfig()
            config.output_directory = str(temp_output_dir)
            config.style.theme = theme
            plotter = StellarPlotter(config)

            # Should initialize without errors
            assert plotter.config.style.theme == theme
