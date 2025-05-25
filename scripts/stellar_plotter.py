#!/usr/bin/env python3
"""
Unified Stellar Structure Plotting System

A comprehensive tool for visualizing stellar structure data from TOV equation solutions.
Supports multiple equation of state (EOS) types and various plot configurations.

Author: Karan Amit Kinariwala
Date: 2025-01-09
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
import math
import logging
import argparse
import yaml
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import sys


# ============================================================================
# CONFIGURATION AND DATA MODELS
# ============================================================================

class PlotType(Enum):
    """Available plot types for stellar structure visualization."""
    SINGLE_PROFILE = "profile"      # Mass + pressure vs radius for single star
    MASS_RADIUS = "mass_radius"     # M-R relation from multiple models
    MASS_DENSITY = "mass_density"   # Mass vs central density
    COMPARATIVE = "comparative"     # Multi-EOS overlay plots


class EOSType(Enum):
    """Supported equation of state types."""
    HYBRID = "hybrid"
    NEUTRON_RELATIVISTIC = "neutron_relativistic"
    ELECTRON_RELATIVISTIC = "electron_relativistic"
    ELECTRON_NON_RELATIVISTIC = "electron_non_relativistic"
    NEUTRON_NON_RELATIVISTIC = "neutron_non_relativistic"


@dataclass
class PlotStyle:
    """Configuration for plot appearance and styling."""
    colors: Dict[str, str] = field(default_factory=lambda: {
        'hybrid': 'blue',
        'neutron_relativistic': 'green', 
        'electron_relativistic': 'red',
        'electron_non_relativistic': 'orange',
        'neutron_non_relativistic': 'purple'
    })
    markers: Dict[str, str] = field(default_factory=lambda: {
        'hybrid': 'o',
        'neutron_relativistic': 's',
        'electron_relativistic': '^',
        'electron_non_relativistic': 'v',
        'neutron_non_relativistic': 'D'
    })
    line_styles: Dict[str, str] = field(default_factory=lambda: {
        'mass': '-',
        'pressure': '--',
        'density': ':'
    })
    figure_size: Tuple[int, int] = (10, 6)
    dpi: int = 300
    grid: bool = True
    legend_location: str = 'best'


@dataclass
class StellarModel:
    """Data model for a single stellar structure solution."""
    mass_solar: float          # Final mass in solar masses
    radius_km: float           # Final radius in kilometers  
    central_density: float     # Central density in g/cm³
    log_central_density: float # log10(central density)
    eos_type: str             # Equation of state identifier
    filename: str             # Source CSV filename
    
    # Optional profile data for single-star plots
    log_r: Optional[np.ndarray] = None    # log10(radius) profile
    log_m: Optional[np.ndarray] = None    # log10(mass) profile  
    log_p: Optional[np.ndarray] = None    # log10(pressure) profile


@dataclass
class EOSDataset:
    """Collection of stellar models for a specific EOS type."""
    eos_type: str
    models: List[StellarModel]
    
    @property
    def masses(self) -> List[float]:
        """Extract masses from all models."""
        return [model.mass_solar for model in self.models]
    
    @property 
    def radii(self) -> List[float]:
        """Extract radii from all models."""
        return [model.radius_km for model in self.models]
    
    @property
    def log_densities(self) -> List[float]:
        """Extract log central densities from all models."""
        return [model.log_central_density for model in self.models]


@dataclass
class PlotConfig:
    """Configuration for plotting operations."""
    data_directory: str = "../data/"
    output_directory: str = "../plots/"
    style: PlotStyle = field(default_factory=PlotStyle)
    
    # EOS-specific file patterns and parsing
    eos_patterns: Dict[str, str] = field(default_factory=lambda: {
        'hybrid': 'tov_solution_magnetic_bps_bbp_polytrope_*.csv',
        'neutron_relativistic': 'neutron_relativistic_rhoc_*.csv',
        'electron_relativistic': 'electron_relativistic_rhoc_*.csv',
        'electron_non_relativistic': 'electron_non_relativistic_rhoc_*.csv',
        'neutron_non_relativistic': 'neutron_non_relativistic_rhoc_*.csv'
    })


# ============================================================================
# CONFIGURATION MANAGER
# ============================================================================

class ConfigManager:
    """Manages configuration loading, validation, and defaults."""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_file: Optional path to YAML configuration file
        """
        self.logger = logging.getLogger(__name__)
        self.config = self._load_config(config_file)
    
    def _load_config(self, config_file: Optional[str]) -> PlotConfig:
        """Load configuration from file or use defaults."""
        if config_file and os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    config_data = yaml.safe_load(f)
                self.logger.info(f"Loaded configuration from {config_file}")
                return self._parse_config_data(config_data)
            except Exception as e:
                self.logger.warning(f"Failed to load config file {config_file}: {e}")
                self.logger.info("Using default configuration")
        
        return PlotConfig()
    
    def _parse_config_data(self, config_data: Dict) -> PlotConfig:
        """Parse configuration data from YAML into PlotConfig object."""
        # This would be expanded to handle full YAML parsing
        # For now, return defaults with any overrides
        config = PlotConfig()
        
        if 'data_directory' in config_data:
            config.data_directory = config_data['data_directory']
        if 'output_directory' in config_data:
            config.output_directory = config_data['output_directory']
            
        return config
    
    def get_config(self) -> PlotConfig:
        """Get the current configuration."""
        return self.config
    
    def validate_paths(self) -> bool:
        """Validate that required directories exist."""
        data_path = Path(self.config.data_directory)
        output_path = Path(self.config.output_directory)
        
        if not data_path.exists():
            self.logger.error(f"Data directory does not exist: {data_path}")
            return False
            
        # Create output directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output directory ready: {output_path}")
        
        return True


# ============================================================================
# EOS DATA PROCESSOR
# ============================================================================

class EOSDataProcessor:
    """Handles CSV file discovery, parsing, and data extraction."""
    
    def __init__(self, config: PlotConfig):
        """
        Initialize the data processor.
        
        Args:
            config: Plot configuration containing paths and patterns
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def discover_files(self, eos_types: List[str]) -> Dict[str, List[str]]:
        """
        Discover CSV files for specified EOS types.
        
        Args:
            eos_types: List of EOS type identifiers
            
        Returns:
            Dictionary mapping EOS type to list of matching file paths
        """
        discovered_files = {}
        
        for eos_type in eos_types:
            if eos_type not in self.config.eos_patterns:
                self.logger.warning(f"Unknown EOS type: {eos_type}")
                continue
                
            pattern = os.path.join(self.config.data_directory, 
                                 self.config.eos_patterns[eos_type])
            files = glob.glob(pattern)
            
            if files:
                discovered_files[eos_type] = sorted(files)
                self.logger.info(f"Found {len(files)} files for {eos_type}")
            else:
                self.logger.warning(f"No files found for {eos_type} with pattern: {pattern}")
        
        return discovered_files
    
    def parse_central_density(self, filename: str, eos_type: str) -> Optional[float]:
        """
        Parse central density from filename based on EOS type.
        
        Args:
            filename: Path to CSV file
            eos_type: EOS type identifier for parsing strategy
            
        Returns:
            Central density in g/cm³, or None if parsing fails
        """
        base = os.path.basename(filename)
        
        try:
            if eos_type == 'hybrid':
                # Pattern: tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv
                match = re.search(r'polytrope_(.*)\.csv', base)
                if match:
                    rho_str = match.group(1)  # e.g. '1.00e+16'
                    return float(rho_str)
            
            else:
                # Pattern: neutron_relativistic_rhoc_5.00pp18.csv
                # Convert encoded notation: '5.00pp18' -> '5.00e+18'
                match = re.search(rf'{eos_type}_rhoc_(.*)\.csv', base)
                if match:
                    rho_str = match.group(1)  # e.g. '5.00pp18'
                    rho_str = rho_str.replace('pp', 'e+').replace('p-', 'e-')
                    # Handle single 'p' for 'e' replacement
                    if 'e' not in rho_str and 'p' in rho_str:
                        rho_str = rho_str.replace('p', 'e', 1)
                    return float(rho_str)
                    
        except (ValueError, AttributeError) as e:
            self.logger.error(f"Failed to parse density from {filename}: {e}")
        
        return None
    
    def load_stellar_model(self, filepath: str, eos_type: str, 
                          load_profile: bool = False) -> Optional[StellarModel]:
        """
        Load stellar model data from CSV file.
        
        Args:
            filepath: Path to CSV file
            eos_type: EOS type identifier
            load_profile: Whether to load full radial profiles
            
        Returns:
            StellarModel object or None if loading fails
        """
        try:
            df = pd.read_csv(filepath)
            
            # Validate required columns
            required_cols = ['log_m[g]', 'log_r[cm]']
            if not all(col in df.columns for col in required_cols):
                self.logger.error(f"Missing required columns in {filepath}")
                return None
            
            # Remove NaN values
            df = df.dropna()
            if df.empty:
                self.logger.error(f"No valid data in {filepath}")
                return None
            
            # Extract final (surface) values
            last_mass_log = df['log_m[g]'].iloc[-1]
            last_radius_log = df['log_r[cm]'].iloc[-1]
            
            # Convert to physical units
            mass_solar = 10**last_mass_log / 1.989e33  # grams to solar masses
            radius_km = 10**last_radius_log / 1e5      # cm to km
            
            # Parse central density
            central_density = self.parse_central_density(filepath, eos_type)
            if central_density is None:
                self.logger.error(f"Could not parse central density from {filepath}")
                return None
            
            log_central_density = math.log10(central_density)
            
            # Create stellar model
            model = StellarModel(
                mass_solar=mass_solar,
                radius_km=radius_km,
                central_density=central_density,
                log_central_density=log_central_density,
                eos_type=eos_type,
                filename=os.path.basename(filepath)
            )
            
            # Load full profiles if requested
            if load_profile:
                model.log_r = df['log_r[cm]'].values
                model.log_m = df['log_m[g]'].values
                if 'log_P[dyne/cm^2]' in df.columns:
                    model.log_p = df['log_P[dyne/cm^2]'].values
            
            return model
            
        except Exception as e:
            self.logger.error(f"Failed to load {filepath}: {e}")
            return None
    
    def load_eos_dataset(self, eos_type: str, load_profiles: bool = False) -> Optional[EOSDataset]:
        """
        Load complete dataset for an EOS type.
        
        Args:
            eos_type: EOS type identifier
            load_profiles: Whether to load full radial profiles
            
        Returns:
            EOSDataset object or None if no valid data found
        """
        files = self.discover_files([eos_type]).get(eos_type, [])
        if not files:
            self.logger.warning(f"No files found for EOS type: {eos_type}")
            return None
        
        models = []
        for filepath in files:
            model = self.load_stellar_model(filepath, eos_type, load_profiles)
            if model:
                models.append(model)
        
        if not models:
            self.logger.error(f"No valid models loaded for EOS type: {eos_type}")
            return None
        
        self.logger.info(f"Loaded {len(models)} models for {eos_type}")
        return EOSDataset(eos_type=eos_type, models=models)


# ============================================================================
# LOGGING SETUP
# ============================================================================

def setup_logging(level: str = "INFO") -> None:
    """Configure logging for the application."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# ============================================================================
# STELLAR PLOTTER ENGINE
# ============================================================================

class StellarPlotter:
    """Handles all plotting functionality with configurable styles."""
    
    def __init__(self, config: PlotConfig):
        """
        Initialize the stellar plotter.
        
        Args:
            config: Plot configuration containing styles and settings
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Set matplotlib defaults
        plt.rcParams['figure.figsize'] = self.config.style.figure_size
        plt.rcParams['figure.dpi'] = self.config.style.dpi
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.grid'] = self.config.style.grid
    
    def plot_single_profile(self, model: StellarModel, output_file: Optional[str] = None) -> str:
        """
        Create mass and pressure profile plot for a single stellar model.
        
        Equivalent to the original plot_m_r_tov.py functionality.
        
        Args:
            model: StellarModel with profile data loaded
            output_file: Optional custom output filename
            
        Returns:
            Path to saved plot file
        """
        if model.log_r is None or model.log_m is None or model.log_p is None:
            raise ValueError("Model must have profile data loaded for single profile plot")
        
        # Convert from logarithmic to linear scale
        r_cm = 10 ** model.log_r
        m_g = 10 ** model.log_m
        p = 10 ** model.log_p
        
        # Unit conversions
        r_km = r_cm / 1e5  # cm to km
        M_sun = 1.989e33   # Solar mass in g
        m_solar = m_g / M_sun  # Mass in solar masses
        
        # Create figure with dual y-axes
        fig, ax1 = plt.subplots(figsize=self.config.style.figure_size)
        ax2 = ax1.twinx()
        
        # Plot mass on primary (left) y-axis
        color_mass = 'red'
        ax1.plot(r_km, m_solar, color=color_mass, linestyle=self.config.style.line_styles['mass'], 
                label='Mass Profile', linewidth=2)
        ax1.set_xlabel("Radius [km]")
        ax1.set_ylabel("Mass [M_sun]", color=color_mass)
        ax1.tick_params(axis='y', colors=color_mass)
        
        # Plot pressure on secondary (right) y-axis
        color_pressure = 'blue'
        ax2.plot(r_km, p, color=color_pressure, linestyle=self.config.style.line_styles['pressure'],
                label='Pressure Profile', linewidth=2)
        ax2.set_ylabel("Pressure [dyne/cm^2]", color=color_pressure)
        ax2.tick_params(axis='y', colors=color_pressure)
        
        # Add title and grid
        plt.title(f"Mass and Pressure Profiles - {model.eos_type.replace('_', ' ').title()}")
        if self.config.style.grid:
            ax1.grid(True, alpha=0.3)
        
        # Create combined legend
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc=self.config.style.legend_location)
        
        plt.tight_layout()
        
        # Generate output filename
        if output_file is None:
            density_str = f"{model.central_density:.2e}".replace('+', 'p').replace('e', 'p')
            output_file = f"profile_{model.eos_type}_rhoc_{density_str}.png"
        
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Saved single profile plot: {output_path}")
        return output_path
    
    def plot_mass_radius_relation(self, datasets: List[EOSDataset], 
                                 output_file: str = "mass_radius_relation.png") -> str:
        """
        Create mass-radius relation plot for multiple EOS datasets.
        
        Equivalent to functionality from plot_mass_radius_relations.py and 
        plot_mass_radius_relations_overplot.py.
        
        Args:
            datasets: List of EOSDataset objects to plot
            output_file: Output filename
            
        Returns:
            Path to saved plot file
        """
        plt.figure(figsize=self.config.style.figure_size)
        
        for dataset in datasets:
            if not dataset.models:
                self.logger.warning(f"No models in dataset {dataset.eos_type}")
                continue
            
            color = self.config.style.colors.get(dataset.eos_type, 'black')
            marker = self.config.style.markers.get(dataset.eos_type, 'o')
            label = dataset.eos_type.replace('_', ' ').title()
            
            plt.scatter(dataset.radii, dataset.masses, 
                       color=color, marker=marker, s=30, alpha=0.7, label=label)
        
        plt.xlabel('Radius (km)')
        plt.ylabel('Mass (Solar Masses)')
        plt.title('Mass-Radius Relations')
        
        if self.config.style.grid:
            plt.grid(True, alpha=0.3)
        
        plt.legend(loc=self.config.style.legend_location)
        plt.tight_layout()
        
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Saved mass-radius relation plot: {output_path}")
        return output_path
    
    def plot_mass_density_relation(self, datasets: List[EOSDataset],
                                  output_file: str = "mass_density_relation.png") -> str:
        """
        Create mass vs central density plot for multiple EOS datasets.
        
        Args:
            datasets: List of EOSDataset objects to plot
            output_file: Output filename
            
        Returns:
            Path to saved plot file
        """
        plt.figure(figsize=self.config.style.figure_size)
        
        for dataset in datasets:
            if not dataset.models:
                self.logger.warning(f"No models in dataset {dataset.eos_type}")
                continue
            
            color = self.config.style.colors.get(dataset.eos_type, 'black')
            marker = self.config.style.markers.get(dataset.eos_type, 'o')
            label = dataset.eos_type.replace('_', ' ').title()
            
            plt.scatter(dataset.log_densities, dataset.masses,
                       color=color, marker=marker, s=30, alpha=0.7, label=label)
        
        plt.xlabel('log10(rho_c) [g/cm^3]')
        plt.ylabel('Mass (Solar Masses)')
        plt.title('Mass vs. log10(Central Density)')
        
        if self.config.style.grid:
            plt.grid(True, alpha=0.3)
        
        plt.legend(loc=self.config.style.legend_location)
        plt.tight_layout()
        
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Saved mass-density relation plot: {output_path}")
        return output_path
    
    def plot_comparative_analysis(self, datasets: List[EOSDataset],
                                 plot_types: List[str] = ["mass_radius", "mass_density"],
                                 output_prefix: str = "comparative") -> List[str]:
        """
        Create comparative plots for multiple EOS types.
        
        Args:
            datasets: List of EOSDataset objects to compare
            plot_types: Types of plots to create ("mass_radius", "mass_density")
            output_prefix: Prefix for output filenames
            
        Returns:
            List of paths to saved plot files
        """
        output_files = []
        
        if "mass_radius" in plot_types:
            output_file = f"{output_prefix}_mass_radius.png"
            path = self.plot_mass_radius_relation(datasets, output_file)
            output_files.append(path)
        
        if "mass_density" in plot_types:
            output_file = f"{output_prefix}_mass_density.png"
            path = self.plot_mass_density_relation(datasets, output_file)
            output_files.append(path)
        
        return output_files
    
    def create_summary_statistics(self, datasets: List[EOSDataset]) -> Dict[str, Dict[str, float]]:
        """
        Calculate summary statistics for each EOS dataset.
        
        Args:
            datasets: List of EOSDataset objects
            
        Returns:
            Dictionary with statistics for each EOS type
        """
        stats = {}
        
        for dataset in datasets:
            if not dataset.models:
                continue
                
            masses = np.array(dataset.masses)
            radii = np.array(dataset.radii)
            densities = np.array(dataset.log_densities)
            
            stats[dataset.eos_type] = {
                'num_models': len(dataset.models),
                'max_mass': float(np.max(masses)),
                'min_mass': float(np.min(masses)),
                'mean_mass': float(np.mean(masses)),
                'max_radius': float(np.max(radii)),
                'min_radius': float(np.min(radii)),
                'mean_radius': float(np.mean(radii)),
                'density_range_log': float(np.max(densities) - np.min(densities))
            }
        
        return stats


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def create_argument_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Unified Stellar Structure Plotting System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot single stellar profile
  python stellar_plotter.py profile --file data/hybrid_1e15.csv --eos-type hybrid
  
  # Create mass-radius relation for hybrid EOS
  python stellar_plotter.py mass-radius --eos-types hybrid
  
  # Compare multiple EOS types
  python stellar_plotter.py compare --eos-types hybrid,neutron_relativistic
  
  # Generate all plots with custom config
  python stellar_plotter.py all --config custom_config.yaml
        """
    )
    
    parser.add_argument('command', choices=['profile', 'mass-radius', 'mass-density', 'compare', 'all'],
                       help='Type of plot to create')
    
    parser.add_argument('--eos-types', type=str, default='hybrid',
                       help='Comma-separated list of EOS types (default: hybrid)')
    
    parser.add_argument('--file', type=str,
                       help='Specific CSV file for single profile plot')
    
    parser.add_argument('--config', type=str, default='stellar_plotter_config.yaml',
                       help='Configuration file path (default: stellar_plotter_config.yaml)')
    
    parser.add_argument('--output-dir', type=str,
                       help='Output directory (overrides config)')
    
    parser.add_argument('--output-prefix', type=str, default='stellar',
                       help='Output filename prefix (default: stellar)')
    
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')
    
    parser.add_argument('--stats', action='store_true',
                       help='Print summary statistics')
    
    return parser


def execute_plotting_command(args) -> None:
    """Execute the plotting command based on parsed arguments."""
    # Setup logging
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)
    
    # Load configuration
    config_manager = ConfigManager(args.config if os.path.exists(args.config) else None)
    config = config_manager.get_config()
    
    # Override output directory if specified
    if args.output_dir:
        config.output_directory = args.output_dir
    
    # Validate configuration
    if not config_manager.validate_paths():
        logger.error("Configuration validation failed")
        return
    
    # Initialize components
    processor = EOSDataProcessor(config)
    plotter = StellarPlotter(config)
    
    # Parse EOS types
    eos_types = [eos.strip() for eos in args.eos_types.split(',')]
    
    try:
        if args.command == 'profile':
            # Single profile plot
            if not args.file:
                logger.error("--file argument required for profile command")
                return
            
            if not os.path.exists(args.file):
                logger.error(f"File not found: {args.file}")
                return
            
            # Determine EOS type from filename or use first specified type
            eos_type = eos_types[0] if eos_types else 'hybrid'
            model = processor.load_stellar_model(args.file, eos_type, load_profile=True)
            
            if model:
                output_file = plotter.plot_single_profile(model)
                logger.info(f"Profile plot saved: {output_file}")
            else:
                logger.error("Failed to load stellar model")
        
        elif args.command in ['mass-radius', 'mass-density', 'compare']:
            # Load datasets for specified EOS types
            datasets = []
            for eos_type in eos_types:
                dataset = processor.load_eos_dataset(eos_type, load_profiles=False)
                if dataset:
                    datasets.append(dataset)
                else:
                    logger.warning(f"No data found for EOS type: {eos_type}")
            
            if not datasets:
                logger.error("No valid datasets loaded")
                return
            
            # Create plots based on command
            if args.command == 'mass-radius':
                output_file = f"{args.output_prefix}_mass_radius.png"
                path = plotter.plot_mass_radius_relation(datasets, output_file)
                logger.info(f"Mass-radius plot saved: {path}")
            
            elif args.command == 'mass-density':
                output_file = f"{args.output_prefix}_mass_density.png"
                path = plotter.plot_mass_density_relation(datasets, output_file)
                logger.info(f"Mass-density plot saved: {path}")
            
            elif args.command == 'compare':
                paths = plotter.plot_comparative_analysis(datasets, 
                                                        ["mass_radius", "mass_density"],
                                                        args.output_prefix)
                for path in paths:
                    logger.info(f"Comparative plot saved: {path}")
            
            # Print statistics if requested
            if args.stats:
                stats = plotter.create_summary_statistics(datasets)
                print("\n" + "="*60)
                print("SUMMARY STATISTICS")
                print("="*60)
                for eos_type, eos_stats in stats.items():
                    print(f"\n{eos_type.replace('_', ' ').title()}:")
                    print(f"  Models: {eos_stats['num_models']}")
                    print(f"  Mass range: {eos_stats['min_mass']:.3f} - {eos_stats['max_mass']:.3f} M☉")
                    print(f"  Radius range: {eos_stats['min_radius']:.1f} - {eos_stats['max_radius']:.1f} km")
                    print(f"  Density range: {eos_stats['density_range_log']:.1f} orders of magnitude")
        
        elif args.command == 'all':
            # Generate all plot types
            datasets = []
            for eos_type in eos_types:
                dataset = processor.load_eos_dataset(eos_type, load_profiles=False)
                if dataset:
                    datasets.append(dataset)
            
            if datasets:
                # Create all comparative plots
                paths = plotter.plot_comparative_analysis(datasets, 
                                                        ["mass_radius", "mass_density"],
                                                        args.output_prefix)
                
                # Print statistics
                stats = plotter.create_summary_statistics(datasets)
                print("\n" + "="*60)
                print("SUMMARY STATISTICS")
                print("="*60)
                for eos_type, eos_stats in stats.items():
                    print(f"\n{eos_type.replace('_', ' ').title()}:")
                    print(f"  Models: {eos_stats['num_models']}")
                    print(f"  Mass range: {eos_stats['min_mass']:.3f} - {eos_stats['max_mass']:.3f} M☉")
                    print(f"  Radius range: {eos_stats['min_radius']:.1f} - {eos_stats['max_radius']:.1f} km")
                
                logger.info(f"Generated {len(paths)} plots successfully")
            else:
                logger.error("No valid datasets found")
    
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        raise


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """Main entry point with CLI support."""
    # Check if command line arguments are provided
    if len(sys.argv) > 1:
        # Use CLI interface
        parser = create_argument_parser()
        args = parser.parse_args()
        execute_plotting_command(args)
    else:
        # Backward compatibility - run Sprint 1 foundation test
        setup_logging()
        logger = logging.getLogger(__name__)
        
        logger.info("Stellar Plotter v2.0 - Sprint 2 Complete")
        logger.info("No command specified - running foundation test")
        
        # Test configuration file loading
        config_file = "stellar_plotter_config.yaml"
        if os.path.exists(config_file):
            logger.info(f"Testing configuration file: {config_file}")
            config_manager = ConfigManager(config_file)
        else:
            logger.info("Using default configuration (no config file found)")
            config_manager = ConfigManager()
        
        config = config_manager.get_config()
        
        if config_manager.validate_paths():
            logger.info("Configuration validation passed")
            
            # Test data processor initialization
            processor = EOSDataProcessor(config)
            discovered = processor.discover_files(['hybrid', 'neutron_relativistic'])
            
            for eos_type, files in discovered.items():
                logger.info(f"EOS {eos_type}: {len(files)} files discovered")
                
            # Test loading a single model to validate the full pipeline
            if 'hybrid' in discovered and discovered['hybrid']:
                logger.info("Testing stellar model loading...")
                first_file = discovered['hybrid'][0]
                model = processor.load_stellar_model(first_file, 'hybrid', load_profile=False)
                if model:
                    logger.info(f"Successfully loaded model: M={model.mass_solar:.3f} M☉, R={model.radius_km:.1f} km, ρc={model.central_density:.2e} g/cm³")
                else:
                    logger.error("Failed to load test model")
            
            logger.info("\nTo use the plotting functionality, run with commands like:")
            logger.info("  python stellar_plotter.py mass-radius --eos-types hybrid")
            logger.info("  python stellar_plotter.py compare --eos-types hybrid,neutron_relativistic")
            logger.info("  python stellar_plotter.py all --stats")
        else:
            logger.error("Configuration validation failed")


if __name__ == "__main__":
    main() 