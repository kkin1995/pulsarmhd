"""
Unit tests for CLI interface.

Tests the command-line interface functionality, argument parsing,
command validation, and help system of the stellar plotter.
"""

import pytest
import sys
import argparse
from unittest.mock import patch, MagicMock
from io import StringIO
from pathlib import Path

# Import the module under test
sys.path.append('../../scripts')
from stellar_plotter import main, create_argument_parser, execute_plotting_command, EOSType


@pytest.mark.unit
class TestCLIArgumentParser:
    """Test suite for CLI argument parser."""

    def test_create_argument_parser(self):
        """Test creation of argument parser."""
        parser = create_argument_parser()

        assert isinstance(parser, argparse.ArgumentParser)
        assert 'Stellar Structure Plotter' in parser.description or 'stellar' in parser.description.lower()

    def test_parser_help_command(self):
        """Test parser help command."""
        parser = create_argument_parser()

        # Capture help output
        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            with pytest.raises(SystemExit):
                parser.parse_args(['--help'])

            help_output = mock_stdout.getvalue()
            assert 'profile' in help_output
            assert 'mass-radius' in help_output

    def test_parser_profile_command(self):
        """Test profile command parsing."""
        parser = create_argument_parser()

        args = parser.parse_args(['profile', '--file', 'test.csv'])

        assert args.command == 'profile'
        assert args.file == 'test.csv'

    def test_parser_mass_radius_command(self):
        """Test mass-radius command parsing."""
        parser = create_argument_parser()

        args = parser.parse_args(['mass-radius', '--eos-types', 'hybrid'])

        assert args.command == 'mass-radius'
        assert args.eos_types == ['hybrid']

    def test_parser_mass_density_command(self):
        """Test mass-density command parsing."""
        parser = create_argument_parser()

        args = parser.parse_args(['mass-density', '--eos-types', 'hybrid', 'neutron_relativistic'])

        assert args.command == 'mass-density'
        assert args.eos_types == ['hybrid', 'neutron_relativistic']

    def test_parser_compare_command(self):
        """Test compare command parsing."""
        parser = create_argument_parser()

        args = parser.parse_args(['compare', '--eos-types', 'hybrid', 'neutron_relativistic'])

        assert args.command == 'compare'
        assert args.eos_types == ['hybrid', 'neutron_relativistic']

    def test_parser_all_command(self):
        """Test all command parsing."""
        parser = create_argument_parser()

        args = parser.parse_args(['all'])

        assert args.command == 'all'

    def test_parser_optional_arguments(self):
        """Test parsing of optional arguments."""
        parser = create_argument_parser()

        args = parser.parse_args([
            'mass-radius',
            '--eos-types', 'hybrid',
            '--theme', 'dark',
            '--dpi', '600',
            '--figure-size', '12', '8',
            '--log-level', 'DEBUG',
            '--stats'
        ])

        assert args.command == 'mass-radius'
        assert args.eos_types == ['hybrid']
        assert args.theme == 'dark'
        assert args.dpi == 600
        assert args.figure_size == [12, 8]
        assert args.log_level == 'DEBUG'
        assert args.stats is True

    def test_parser_invalid_command(self):
        """Test parsing with invalid command."""
        parser = create_argument_parser()

        with patch('sys.stderr', new_callable=StringIO):
            with pytest.raises(SystemExit):
                parser.parse_args(['invalid_command'])

    def test_parser_profile_requires_file(self):
        """Test that profile command requires --file argument."""
        parser = create_argument_parser()

        with patch('sys.stderr', new_callable=StringIO):
            with pytest.raises(SystemExit):
                parser.parse_args(['profile'])  # Missing --file

    def test_parser_figure_size_validation(self):
        """Test figure size argument validation."""
        parser = create_argument_parser()

        # Valid figure size
        args = parser.parse_args(['mass-radius', '--figure-size', '10', '6'])
        assert args.figure_size == [10, 6]


@pytest.mark.unit
class TestCLIMainFunction:
    """Test suite for CLI main function."""

    @patch('stellar_plotter.execute_plotting_command')
    def test_main_with_arguments(self, mock_execute):
        """Test main function with command line arguments."""
        test_args = ['stellar_plotter.py', 'mass-radius', '--eos-types', 'hybrid']
        with patch('sys.argv', test_args):
            main()

        # Should call execute_plotting_command
        mock_execute.assert_called_once()

    @patch('stellar_plotter.setup_logging')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.EOSDataProcessor')
    def test_main_without_arguments(self, mock_processor, mock_config_manager, mock_setup_logging):
        """Test main function without command line arguments (foundation test)."""
        # Mock the configuration manager
        mock_config = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock the processor
        mock_processor_instance = MagicMock()
        mock_processor_instance.discover_files.return_value = {'hybrid': ['test_file.csv']}

        # Create a proper mock stellar model with numeric attributes
        mock_model = MagicMock()
        mock_model.mass_solar = 2.0
        mock_model.radius_km = 12.0
        mock_model.central_density = 1.0e16
        mock_processor_instance.load_stellar_model.return_value = mock_model
        mock_processor.return_value = mock_processor_instance

        test_args = ['stellar_plotter.py']  # No additional arguments
        with patch('sys.argv', test_args):
            main()

        # Should run foundation test
        mock_setup_logging.assert_called()
        mock_config_manager.assert_called()


@pytest.mark.unit
class TestCLIExecutePlottingCommand:
    """Test suite for execute_plotting_command function."""

    @patch('stellar_plotter.StellarPlotter')
    @patch('stellar_plotter.EOSDataProcessor')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_execute_profile_command(self, mock_setup_logging, mock_config_manager, mock_processor, mock_plotter):
        """Test executing profile command."""
        # Mock configuration
        mock_config = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock processor with proper stellar model attributes
        mock_model = MagicMock()
        mock_model.eos_type = "hybrid"
        mock_model.central_density = 1.0e16
        mock_processor_instance = MagicMock()
        mock_processor_instance.load_stellar_model.return_value = mock_model
        mock_processor.return_value = mock_processor_instance

        # Mock plotter
        mock_plotter_instance = MagicMock()
        mock_plotter_instance.plot_single_profile.return_value = "test_output.png"
        mock_plotter.return_value = mock_plotter_instance

        # Create mock args with proper attributes
        args = MagicMock()
        args.command = 'profile'
        args.file = 'test.csv'
        args.log_level = 'INFO'
        args.output_prefix = ''

        # Mock os.path.exists to return True
        with patch('os.path.exists', return_value=True):
            execute_plotting_command(args)

        # Verify calls
        mock_setup_logging.assert_called_with('INFO')
        mock_processor_instance.load_stellar_model.assert_called()
        mock_plotter_instance.plot_single_profile.assert_called()

    @patch('stellar_plotter.StellarPlotter')
    @patch('stellar_plotter.EOSDataProcessor')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_execute_mass_radius_command(self, mock_setup_logging, mock_config_manager, mock_processor, mock_plotter):
        """Test executing mass-radius command."""
        # Mock configuration
        mock_config = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock processor
        mock_dataset = MagicMock()
        mock_processor_instance = MagicMock()
        mock_processor_instance.load_eos_dataset.return_value = mock_dataset
        mock_processor.return_value = mock_processor_instance

        # Mock plotter
        mock_plotter_instance = MagicMock()
        mock_plotter_instance.plot_mass_radius_relation.return_value = "test_output.png"
        mock_plotter.return_value = mock_plotter_instance

        # Create mock args
        args = MagicMock()
        args.command = 'mass-radius'
        args.eos_types = ['hybrid']
        args.log_level = 'INFO'
        args.output_prefix = ''

        execute_plotting_command(args)

        # Verify calls
        mock_setup_logging.assert_called_with('INFO')
        mock_processor_instance.load_eos_dataset.assert_called_with('hybrid')
        mock_plotter_instance.plot_mass_radius_relation.assert_called()

    @patch('stellar_plotter.StellarPlotter')
    @patch('stellar_plotter.EOSDataProcessor')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_execute_compare_command(self, mock_setup_logging, mock_config_manager, mock_processor, mock_plotter):
        """Test executing compare command."""
        # Mock configuration
        mock_config = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock processor
        mock_dataset = MagicMock()
        mock_processor_instance = MagicMock()
        mock_processor_instance.load_eos_dataset.return_value = mock_dataset
        mock_processor.return_value = mock_processor_instance

        # Mock plotter
        mock_plotter_instance = MagicMock()
        mock_plotter_instance.plot_comparative_analysis.return_value = ["test_output1.png", "test_output2.png"]
        mock_plotter.return_value = mock_plotter_instance

        # Create mock args
        args = MagicMock()
        args.command = 'compare'
        args.eos_types = ['hybrid', 'neutron_relativistic']
        args.log_level = 'INFO'
        args.output_prefix = ''

        execute_plotting_command(args)

        # Verify calls
        mock_setup_logging.assert_called_with('INFO')
        mock_plotter_instance.plot_comparative_analysis.assert_called()

    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_execute_invalid_config_path(self, mock_setup_logging, mock_config_manager):
        """Test executing command with invalid configuration path."""
        mock_config_manager.return_value.validate_paths.return_value = False

        args = MagicMock()
        args.command = 'mass-radius'
        args.log_level = 'INFO'

        # Should handle gracefully and not raise exception
        execute_plotting_command(args)

        mock_setup_logging.assert_called_with('INFO')


@pytest.mark.unit
class TestCLIConfigurationOverrides:
    """Test suite for CLI configuration overrides."""

    @patch('stellar_plotter.StellarPlotter')
    @patch('stellar_plotter.EOSDataProcessor')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_config_override_theme(self, mock_setup_logging, mock_config_manager, mock_processor, mock_plotter):
        """Test theme override from CLI."""
        # Mock configuration
        mock_config = MagicMock()
        mock_config.style = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock processor and plotter
        mock_processor.return_value = MagicMock()
        mock_plotter.return_value = MagicMock()

        # Create mock args with theme override
        args = MagicMock()
        args.command = 'mass-radius'
        args.eos_types = ['hybrid']
        args.theme = 'dark'
        args.log_level = 'INFO'
        args.output_prefix = ''

        execute_plotting_command(args)

        # Config should be updated with CLI override
        assert mock_config.style.theme == 'dark'

    @patch('stellar_plotter.StellarPlotter')
    @patch('stellar_plotter.EOSDataProcessor')
    @patch('stellar_plotter.ConfigManager')
    @patch('stellar_plotter.setup_logging')
    def test_config_override_dpi(self, mock_setup_logging, mock_config_manager, mock_processor, mock_plotter):
        """Test DPI override from CLI."""
        # Mock configuration
        mock_config = MagicMock()
        mock_config.style = MagicMock()
        mock_config_manager.return_value.get_config.return_value = mock_config
        mock_config_manager.return_value.validate_paths.return_value = True

        # Mock processor and plotter
        mock_processor.return_value = MagicMock()
        mock_plotter.return_value = MagicMock()

        # Create mock args with DPI override
        args = MagicMock()
        args.command = 'mass-radius'
        args.eos_types = ['hybrid']
        args.dpi = 600
        args.log_level = 'INFO'
        args.output_prefix = ''

        execute_plotting_command(args)

        # Config should be updated with CLI override
        assert mock_config.style.dpi == 600


@pytest.mark.unit
class TestCLIEOSTypeHandling:
    """Test suite for CLI EOS type handling."""

    def test_eos_type_conversion_single(self):
        """Test conversion of single EOS type string."""
        parser = create_argument_parser()
        args = parser.parse_args(['mass-radius', '--eos-types', 'hybrid'])

        # The CLI should handle string EOS types
        assert 'hybrid' in args.eos_types

    def test_eos_type_conversion_multiple(self):
        """Test conversion of multiple EOS type strings."""
        parser = create_argument_parser()
        args = parser.parse_args(['compare', '--eos-types', 'hybrid', 'neutron_relativistic'])

        assert 'hybrid' in args.eos_types
        assert 'neutron_relativistic' in args.eos_types

    def test_eos_type_validation(self):
        """Test validation of EOS type strings."""
        parser = create_argument_parser()

        # Valid EOS types
        valid_types = ['hybrid', 'neutron_relativistic', 'electron_relativistic',
                      'electron_non_relativistic', 'neutron_non_relativistic']

        for eos_type in valid_types:
            args = parser.parse_args(['mass-radius', '--eos-types', eos_type])
            assert eos_type in args.eos_types


@pytest.mark.unit
class TestCLILogging:
    """Test suite for CLI logging configuration."""

    @patch('stellar_plotter.setup_logging')
    def test_logging_level_debug(self, mock_setup_logging):
        """Test DEBUG logging level configuration."""
        args = MagicMock()
        args.command = 'mass-radius'
        args.log_level = 'DEBUG'

        # Mock other dependencies
        with patch('stellar_plotter.ConfigManager') as mock_config_manager:
            mock_config_manager.return_value.validate_paths.return_value = False

            execute_plotting_command(args)

        # Logging should be configured with DEBUG level
        mock_setup_logging.assert_called_with('DEBUG')

    @patch('stellar_plotter.setup_logging')
    def test_logging_level_info(self, mock_setup_logging):
        """Test INFO logging level configuration."""
        args = MagicMock()
        args.command = 'mass-radius'
        args.log_level = 'INFO'

        # Mock other dependencies
        with patch('stellar_plotter.ConfigManager') as mock_config_manager:
            mock_config_manager.return_value.validate_paths.return_value = False

            execute_plotting_command(args)

        # Logging should be configured with INFO level
        mock_setup_logging.assert_called_with('INFO')
