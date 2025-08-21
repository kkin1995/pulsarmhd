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
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import glob
import os
import re
import math
import logging
import argparse
import yaml
from typing import Dict, List, Optional, Tuple
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
    figure_size: Tuple[float, float] = (10.0, 6.0)
    dpi: int = 300
    grid: bool = True
    legend_location: str = 'best'
    theme: str = 'publication'  # 'publication', 'presentation', 'dark', 'colorblind'
    font_family: str = 'serif'  # 'serif', 'sans-serif', 'monospace'
    line_width: float = 2.0
    marker_size: float = 30.0
    alpha: float = 0.7
    grid_alpha: float = 0.3
    show_title: bool = True
    show_legend: bool = True
    connect_points: bool = False       # draw a line through sorted points
    curve_width: float = 1.2
    curve_alpha: float = 1.0
    mark_max_point: bool = False


class PlotTheme:
    """Predefined plot themes for different use cases."""
    
    @staticmethod
    def get_theme_config(theme_name: str) -> Dict[str, any]:
        """Get matplotlib configuration for a specific theme."""
        themes = {
            'publication': {
                'font.family': 'serif',
                'font.size': 12,
                'axes.linewidth': 1.2,
                'axes.spines.top': False,
                'axes.spines.right': False,
                'axes.grid': True,
                'grid.alpha': 0.3,
                'legend.frameon': False,
                'legend.fancybox': False,
                'legend.shadow': False,
                'figure.facecolor': 'white',
                'axes.facecolor': 'white'
            },
            'presentation': {
                'font.family': 'sans-serif',
                'font.size': 14,
                'axes.linewidth': 2.0,
                'axes.spines.top': False,
                'axes.spines.right': False,
                'axes.grid': True,
                'grid.alpha': 0.4,
                'legend.frameon': True,
                'legend.fancybox': True,
                'figure.facecolor': 'white',
                'axes.facecolor': '#f8f8f8'
            },
            'dark': {
                'font.family': 'sans-serif',
                'font.size': 12,
                'text.color': 'white',
                'axes.labelcolor': 'white',
                'axes.edgecolor': 'white',
                'xtick.color': 'white',
                'ytick.color': 'white',
                'axes.spines.top': False,
                'axes.spines.right': False,
                'axes.grid': True,
                'grid.color': 'gray',
                'grid.alpha': 0.3,
                'legend.frameon': True,
                'legend.facecolor': '#2e2e2e',
                'legend.edgecolor': 'white',
                'figure.facecolor': '#1e1e1e',
                'axes.facecolor': '#2e2e2e'
            },
            'colorblind': {
                'font.family': 'sans-serif',
                'font.size': 12,
                'axes.linewidth': 1.2,
                'axes.spines.top': False,
                'axes.spines.right': False,
                'axes.grid': True,
                'grid.alpha': 0.3,
                'legend.frameon': True,
                'figure.facecolor': 'white',
                'axes.facecolor': 'white'
            }
        }
        
        return themes.get(theme_name, themes['publication'])
    
    @staticmethod
    def get_colorblind_palette() -> Dict[str, str]:
        """Get colorblind-friendly color palette."""
        return {
            'hybrid': '#1f77b4',           # Blue
            'neutron_relativistic': '#ff7f0e',     # Orange  
            'electron_relativistic': '#2ca02c',    # Green
            'electron_non_relativistic': '#d62728', # Red
            'neutron_non_relativistic': '#9467bd'   # Purple
        }


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
class DataValidation:
    """Data validation and quality checks."""
    min_mass_solar: float = 0.01      # Minimum physical mass (solar masses)
    max_mass_solar: float = 10.0      # Maximum reasonable mass (solar masses)
    min_radius_km: float = 1.0        # Minimum physical radius (km)
    max_radius_km: float = 2000.0     # Maximum reasonable radius (km) - increased for relativistic polytrope
    min_density_log: float = 10.0     # Minimum log density (g/cm³)
    max_density_log: float = 20.0     # Maximum log density (g/cm³)
    skip_r_end_artifacts: bool = True       # enable by default
    r_end_km: float = 1000.0                # C++ uses r_end = 1e8 cm => 1000 km
    r_end_tolerance_km: float = 1.0         # how close to 1000 km counts as "ended at r_end"
    
    def validate_model(self, model: StellarModel) -> Tuple[bool, List[str]]:
        """
        Validate a stellar model for physical reasonableness.
        
        Args:
            model: StellarModel to validate
            
        Returns:
            Tuple of (is_valid, list_of_warnings)
        """
        warnings = []
        is_valid = True
        
        # Mass validation
        if model.mass_solar < self.min_mass_solar:
            warnings.append(f"Mass {model.mass_solar:.3f} M☉ below minimum {self.min_mass_solar}")
            is_valid = False
        elif model.mass_solar > self.max_mass_solar:
            warnings.append(f"Mass {model.mass_solar:.3f} M☉ above maximum {self.max_mass_solar}")
            is_valid = False
        
        # Radius validation
        if model.radius_km < self.min_radius_km:
            warnings.append(f"Radius {model.radius_km:.1f} km below minimum {self.min_radius_km}")
            is_valid = False
        elif model.radius_km > self.max_radius_km:
            warnings.append(f"Radius {model.radius_km:.1f} km above maximum {self.max_radius_km}")
            is_valid = False
        
        # Density validation
        if model.log_central_density < self.min_density_log:
            warnings.append(f"Log density {model.log_central_density:.1f} below minimum {self.min_density_log}")
            is_valid = False
        elif model.log_central_density > self.max_density_log:
            warnings.append(f"Log density {model.log_central_density:.1f} above maximum {self.max_density_log}")
            is_valid = False
        
        # Physical consistency checks
        if model.eos_type in ['electron_relativistic', 'electron_non_relativistic']:
            # White dwarf mass should be < 1.4 M☉ (Chandrasekhar limit)
            if model.mass_solar > 1.4:
                warnings.append(f"White dwarf mass {model.mass_solar:.3f} M☉ exceeds Chandrasekhar limit")
        
        elif model.eos_type in ['neutron_relativistic', 'neutron_non_relativistic']:
            # Neutron star should have reasonable mass and radius
            if model.mass_solar < 0.1:
                warnings.append(f"Neutron star mass {model.mass_solar:.3f} M☉ unusually low")
            if model.radius_km > 20.0:
                warnings.append(f"Neutron star radius {model.radius_km:.1f} km unusually large")
        
        return is_valid, warnings


@dataclass
class PerformanceConfig:
    """Configuration for performance optimizations."""
    enable_caching: bool = True
    max_cache_size: int = 100
    enable_parallel_loading: bool = True
    chunk_size: int = 10
    memory_limit_mb: int = 1000


@dataclass
class PlotConfig:
    """Configuration for plotting operations."""
    data_directory: str = "./data/"
    output_directory: str = "./plots/"
    style: PlotStyle = field(default_factory=PlotStyle)
    validation: DataValidation = field(default_factory=DataValidation)
    performance: PerformanceConfig = field(default_factory=PerformanceConfig)
    
    # EOS-specific file patterns and parsing
    eos_patterns: Dict[str, str] = field(default_factory=lambda: {
        'hybrid': 'tov_solution_magnetic_bps_bbp_polytrope_*.csv',
        'neutron_relativistic': 'neutron_relativistic_rhoc_*.csv',
        'electron_relativistic': 'electron_relativistic_rhoc_*.csv',
        'electron_non_relativistic': 'electron_non_relativistic_rhoc_*.csv',
        'neutron_non_relativistic': 'neutron_non_relativistic_rhoc_*.csv'
    })
    
    # Enhanced filename parsing with C++ compatibility
    def parse_cpp_encoded_density(self, density_str: str) -> float:
        """
        Parse density string using C++ encoding rules.
        
        C++ get_filename() does:
        1. std::replace('+', 'p') 
        2. std::replace('e', 'p')
        
        So '1.00e+16' becomes '1.00pp16'
        """
        # Reverse the C++ encoding
        if 'pp' in density_str:
            # Handle positive exponent: '1.00pp16' -> '1.00e+16'
            density_str = density_str.replace('pp', 'e+')
        elif density_str.count('p') == 1:
            # Handle negative exponent or single 'e': '1.50p-03' -> '1.50e-03'
            density_str = density_str.replace('p', 'e', 1)
        
        return float(density_str)


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

        # Validation block
        if 'validation' in config_data and config_data['validation'] is not None:
            v = config_data['validation']
            for key in ['min_mass_solar','max_mass_solar','min_radius_km','max_radius_km',
                        'min_density_log','max_density_log',
                        'skip_r_end_artifacts','r_end_km','r_end_tolerance_km']:
                if key in v:
                    setattr(config.validation, key, v[key])

        # Style block (optional but handy)
        if 'style' in config_data and config_data['style'] is not None:
            s = config_data['style']
            if 'theme' in s:               config.style.theme = s['theme']
            if 'figure_size' in s:         config.style.figure_size = tuple(s['figure_size'])
            if 'dpi' in s:                 config.style.dpi = s['dpi']
            if 'grid' in s:                config.style.grid = bool(s['grid'])
            if 'grid_alpha' in s:          config.style.grid_alpha = float(s['grid_alpha'])
            if 'legend_location' in s:     config.style.legend_location = s['legend_location']
            if 'font_family' in s:         config.style.font_family = s['font_family']
            if 'line_width' in s:          config.style.line_width = float(s['line_width'])
            if 'marker_size' in s:         config.style.marker_size = float(s['marker_size'])
            if 'alpha' in s:               config.style.alpha = float(s['alpha'])
            if 'colors' in s and isinstance(s['colors'], dict):
                config.style.colors.update(s['colors'])
            if 'markers' in s and isinstance(s['markers'], dict):
                config.style.markers.update(s['markers'])
            if 'line_styles' in s and isinstance(s['line_styles'], dict):
                config.style.line_styles.update(s['line_styles'])
            if 'show_title' in s:        config.style.show_title = bool(s['show_title'])
            if 'show_legend' in s:       config.style.show_legend = bool(s['show_legend'])
            if 'connect_points' in s:    config.style.connect_points = bool(s['connect_points'])
            if 'curve_width' in s:       config.style.curve_width = float(s['curve_width'])
            if 'curve_alpha' in s:       config.style.curve_alpha = float(s['curve_alpha'])
            if 'mark_max_point' in s:    config.style.mark_max_point = bool(s['mark_max_point'])


        # EOS file patterns (optional override)
        if 'eos_patterns' in config_data and isinstance(config_data['eos_patterns'], dict):
            config.eos_patterns.update(config_data['eos_patterns'])

        # Performance (optional)
        if 'performance' in config_data and config_data['performance'] is not None:
            p = config_data['performance']
            if 'enable_parallel_loading' in p:
                config.performance.enable_parallel_loading = bool(p['enable_parallel_loading'])
            if 'chunk_size' in p:
                config.performance.chunk_size = int(p['chunk_size'])
            
        return config
    
    def get_config(self) -> PlotConfig:
        """Get the current configuration."""
        return self.config
    
    def validate_paths(self) -> bool:
        """Validate that required directories exist."""
        data_path = Path(self.config.data_directory).expanduser().resolve()
        output_path = Path(self.config.output_directory).expanduser().resolve()

        if not data_path.exists():
            self.logger.error(f"Data directory does not exist: {data_path}")
            return False

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
        Uses enhanced C++ compatibility parsing.
        
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
                # Use enhanced C++ compatibility parsing
                match = re.search(rf'{eos_type}_rhoc_(.*)\.csv', base)
                if match:
                    rho_str = match.group(1)  # e.g. '5.00pp18'
                    return self.config.parse_cpp_encoded_density(rho_str)
                    
        except (ValueError, AttributeError) as e:
            self.logger.error(f"Failed to parse density from {filename}: {e}")
        
        return None
    
    def load_stellar_model(self, filepath: str, eos_type: str, 
                          load_profile: bool = False) -> Optional[StellarModel]:
        """
        Load stellar model data from CSV file with enhanced validation.
        
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
            
            # Remove NaN values and validate data quality
            initial_rows = len(df)
            df = df.dropna()
            if df.empty:
                self.logger.error(f"No valid data in {filepath}")
                return None
            
            if len(df) < initial_rows * 0.5:  # More than 50% NaN values
                self.logger.warning(f"High NaN ratio in {filepath}: {initial_rows - len(df)}/{initial_rows} rows removed")
            
            # Extract final (surface) values
            last_mass_log = df['log_m[g]'].iloc[-1]
            last_radius_log = df['log_r[cm]'].iloc[-1]
            
            # Convert to physical units
            mass_solar = 10**last_mass_log / 1.989e33  # grams to solar masses
            radius_km = 10**last_radius_log / 1e5      # cm to km
            
            if self.config.validation and getattr(self.config.validation, 'skip_r_end_artifacts', False):
                r_end_km = getattr(self.config.validation, 'r_end_km', 1000.0)
                tol_km  = getattr(self.config.validation, 'r_end_tolerance_km', 1.0)

                # quick geometric artifact test
                ended_at_r_end = abs(radius_km - r_end_km) <= tol_km

                # optional: strengthen the decision if pressure column is present
                if 'log_P[dyne/cm^2]' in df.columns:
                    last_lp = float(df['log_P[dyne/cm^2]'].iloc[-1])
                    min_lp  = float(df['log_P[dyne/cm^2]'].min())
                    # if we ended at the very last step without ever reaching the true surface,
                    # the final pressure typically isn't the minimum attainable by the EOS trace
                    # (small tolerance in log10 space)
                    not_at_surface = not np.isclose(last_lp, min_lp, atol=5e-4)
                    ended_at_r_end = ended_at_r_end or not_at_surface and (radius_km >= r_end_km - tol_km)

                if ended_at_r_end:
                    self.logger.warning(f"{os.path.basename(filepath)} ended at r_end≈{r_end_km:.0f} km; skipping.")
                    return None

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
            
            # Validate model if validation is enabled
            if self.config.validation:
                is_valid, warnings = self.config.validation.validate_model(model)
                if warnings:
                    for warning in warnings:
                        self.logger.warning(f"{filepath}: {warning}")
                if not is_valid and len(warnings) > 1:  # Allow single warnings
                    self.logger.error(f"Model validation failed for {filepath}")
                    return None
            
            # Load full profiles if requested
            if load_profile:
                model.log_r = df['log_r[cm]'].values
                model.log_m = df['log_m[g]'].values
                if 'log_P[dyne/cm^2]' in df.columns:
                    model.log_p = df['log_P[dyne/cm^2]'].values
                else:
                    self.logger.warning(f"No pressure data in {filepath}")
            
            return model
            
        except Exception as e:
            self.logger.error(f"Failed to load {filepath}: {e}")
            return None
    
    def load_eos_dataset(self, eos_type: str, load_profiles: bool = False) -> Optional[EOSDataset]:
        """
        Load complete dataset for an EOS type with performance optimizations.
        
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
        failed_count = 0
        
        # Performance optimization: process in chunks if enabled
        if self.config.performance.enable_parallel_loading and len(files) > self.config.performance.chunk_size:
            self.logger.info(f"Processing {len(files)} files in chunks of {self.config.performance.chunk_size}")
        
        for i, filepath in enumerate(files):
            model = self.load_stellar_model(filepath, eos_type, load_profiles)
            if model:
                models.append(model)
            else:
                failed_count += 1
            
            # Progress reporting for large datasets
            if len(files) > 20 and (i + 1) % 10 == 0:
                self.logger.info(f"Processed {i + 1}/{len(files)} files for {eos_type}")
        
        if not models:
            self.logger.error(f"No valid models loaded for EOS type: {eos_type}")
            return None
        
        if failed_count > 0:
            self.logger.warning(f"Failed to load {failed_count}/{len(files)} files for {eos_type}")
        
        self.logger.info(f"Successfully loaded {len(models)} models for {eos_type}")
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
    """Handles all plotting functionality with configurable styles and themes."""
    
    def __init__(self, config: PlotConfig):
        """
        Initialize the stellar plotter with enhanced theming.
        
        Args:
            config: Plot configuration containing styles and settings
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Apply theme configuration
        self._apply_theme()
        
        # Set matplotlib defaults
        plt.rcParams['figure.figsize'] = self.config.style.figure_size
        plt.rcParams['figure.dpi'] = self.config.style.dpi
        
        # Apply colorblind palette if requested
        if self.config.style.theme == 'colorblind':
            self.config.style.colors = PlotTheme.get_colorblind_palette()
    
    def _apply_theme(self):
        """Apply the selected theme to matplotlib."""
        theme_config = PlotTheme.get_theme_config(self.config.style.theme)
        plt.rcParams.update(theme_config)
        self.logger.info(f"Applied {self.config.style.theme} theme")
    
    def _get_plot_style(self, eos_type: str) -> Dict[str, any]:
        """Get plotting style parameters for an EOS type."""
        return {
            'color': self.config.style.colors.get(eos_type, 'black'),
            'marker': self.config.style.markers.get(eos_type, 'o'),
            's': self.config.style.marker_size,
            'alpha': self.config.style.alpha,
            'linewidth': self.config.style.line_width
        }
    
    def plot_single_profile(self, model: StellarModel, output_file: Optional[str] = None) -> str:
        """
        Create mass and pressure profile plot for a single stellar model.
        Enhanced with theme support and improved styling.
        
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
        
        # Enhanced styling based on theme
        mass_color = '#d62728' if self.config.style.theme != 'dark' else '#ff6b6b'
        pressure_color = '#1f77b4' if self.config.style.theme != 'dark' else '#4dabf7'
        
        # Plot mass on primary (left) y-axis
        ax1.plot(r_km, m_solar, color=mass_color, 
                linestyle=self.config.style.line_styles['mass'], 
                label='Mass Profile', linewidth=self.config.style.line_width)
        ax1.set_xlabel("Radius [km]")
        ax1.set_ylabel("Mass [M☉]", color=mass_color)
        ax1.tick_params(axis='y', colors=mass_color)
        
        # Plot pressure on secondary (right) y-axis
        ax2.plot(r_km, p, color=pressure_color, 
                linestyle=self.config.style.line_styles['pressure'],
                label='Pressure Profile', linewidth=self.config.style.line_width)
        ax2.set_ylabel("Pressure [dyne/cm²]", color=pressure_color)
        ax2.tick_params(axis='y', colors=pressure_color)
        
        # Enhanced title with model information
        title = f"Stellar Structure Profile - {model.eos_type.replace('_', ' ').title()}\n"
        title += f"M = {model.mass_solar:.3f} M☉, R = {model.radius_km:.1f} km, ρc = {model.central_density:.2e} g/cm³"
        plt.title(title, fontsize=12, pad=20)
        
        # Enhanced grid
        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)
        else:
            ax.grid(False)

        # Create combined legend with enhanced styling
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, 
                  loc=self.config.style.legend_location,
                  framealpha=0.9, fancybox=True, shadow=True)
        
        plt.tight_layout()
        
        # Generate output filename with enhanced encoding
        if output_file is None:
            # Use C++ compatible encoding for consistency
            density_str = f"{model.central_density:.2e}".replace('+', 'p').replace('e', 'p')
            output_file = f"profile_{model.eos_type}_rhoc_{density_str}.png"
        
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight', 
                   facecolor='auto', edgecolor='none')
        plt.close()
        
        self.logger.info(f"Saved single profile plot: {output_path}")
        return output_path
    
    def plot_mass_radius_relation(self, datasets: List[EOSDataset], 
                                 output_file: str = "mass_radius_relation.png") -> str:
        """
        Create enhanced mass-radius relation plot for multiple EOS datasets.
        
        Args:
            datasets: List of EOSDataset objects to plot
            output_file: Output filename
            
        Returns:
            Path to saved plot file
        """
        fig, ax = plt.subplots(figsize=self.config.style.figure_size)
        
        # Always define a cap (∞ if no validation/max set)
        if getattr(self.config, "validation", None) and hasattr(self.config.validation, "max_radius_km"):
            radius_cap = float(self.config.validation.max_radius_km)
        else:
            radius_cap = float("inf")

        for dataset in datasets:
            if not dataset.models:
                self.logger.warning(f"No models in dataset {dataset.eos_type}")
                continue
            
            # filter out unphysical/undesired large-radius points (e.g., 800+ km)
            models_kept = [m for m in dataset.models if m.radius_km <= radius_cap]
            dropped = len(dataset.models) - len(models_kept)
            if dropped > 0:
                self.logger.info(f"{dataset.eos_type}: dropped {dropped} models with R > {radius_cap} km")

            if not models_kept:
                self.logger.warning(f"No models to plot after radius filter for {dataset.eos_type}")
                continue

            # sort by central density for a clean curve
            models_sorted = sorted(models_kept, key=lambda m: m.log_central_density)
            radii  = [m.radius_km  for m in models_sorted]
            masses = [m.mass_solar for m in models_sorted]

            style = self._get_plot_style(dataset.eos_type)
            ax.scatter(
                radii, masses,
                marker=style['marker'], s=style['s'],
                alpha=style['alpha'],
                facecolor=style['color'], edgecolor=style['color'], linewidths=0.0,
                label=dataset.eos_type.replace('_', ' ').title(),
                zorder=2
            )

            if getattr(self.config.style, "connect_points", False):
                ax.plot(
                    radii, masses, "-",
                    lw=self.config.style.curve_width,
                    alpha=self.config.style.curve_alpha,
                    color=style['color'],
                    zorder=1
                )

            if getattr(self.config.style, "mark_max_point", False):
                i_max = int(np.argmax(masses))
                ax.scatter([radii[i_max]], [masses[i_max]], s=50, marker='*',
                        color='0.15', zorder=3)
                ax.annotate(r"$M_{\max}$", xy=(radii[i_max], masses[i_max]),
                            xytext=(6, 6), textcoords="offset points", fontsize=8)
        
        # Labels
        ax.set_xlabel(r"$R\;(\mathrm{km})$")
        ax.set_ylabel(r"$M\;(M_\odot)$")
        
        # Axes cosmetics
        ax.tick_params(direction='in', which='both', top=True, right=True, length=4)
        ax.tick_params(which='minor', length=2)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)

        # Title/legend logic
        if self.config.style.show_title:
            ax.set_title('Mass–Radius Relations for Compact Objects', pad=10)
        if self.config.style.show_legend and len(datasets) > 1:
            ax.legend(loc=self.config.style.legend_location)

        plt.tight_layout()

        output_path = os.path.join(self.config.output_directory, output_file)
        fig.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')
        self.logger.info(f"Saved mass-radius relation plot: {output_path}")

        pdf_path = os.path.splitext(output_path)[0] + ".pdf"
        fig.savefig(pdf_path, bbox_inches='tight')

        plt.close(fig)

        return output_path
    
    def plot_mass_density_relation(self, datasets: List[EOSDataset],
                                  output_file: str = "mass_density_relation.png") -> str:
        """
        Create enhanced mass vs central density plot for multiple EOS datasets.
        
        Args:
            datasets: List of EOSDataset objects to plot
            output_file: Output filename
            
        Returns:
            Path to saved plot file
        """
        fig, ax = plt.subplots(figsize=self.config.style.figure_size)
        
        for dataset in datasets:
            if not dataset.models:
                self.logger.warning(f"No models in dataset {dataset.eos_type}")
                continue
            
            style = self._get_plot_style(dataset.eos_type)
            label = dataset.eos_type.replace('_', ' ').title()
            
            ax.scatter(dataset.log_densities, dataset.masses,
                    marker=style['marker'], s=style['s'],
                    alpha=style['alpha'], facecolor=style['color'],
                    edgecolor=style['color'], linewidths=0.0,
                    label=dataset.eos_type.replace('_',' ').title(), zorder=2)
                    
            if getattr(self.config.style, "connect_points", False):
                xy = sorted(zip(dataset.log_densities, dataset.masses))
                ax.plot([x for x,_ in xy], [y for _,y in xy], "-",
                        lw=self.config.style.curve_width,
                        alpha=self.config.style.curve_alpha,
                        color=style['color'], zorder=1)
        
        ax.set_xlabel(r"$\log_{10}\rho_c\;[\mathrm{g\,cm^{-3}}]$")
        ax.set_ylabel(r"$M\;(M_\odot)$")
        ax.tick_params(direction='in', which='both', top=True, right=True, length=4)
        ax.tick_params(which='minor', length=2)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)
        if self.config.style.show_title:
            ax.set_title('Mass vs. Central Density', pad=10)
        if self.config.style.show_legend and len(datasets) > 1:
            ax.legend(loc=self.config.style.legend_location)
        plt.tight_layout()
        
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight',
                   facecolor='auto', edgecolor='none')
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
    """Create enhanced command-line argument parser with Sprint 3 features."""
    parser = argparse.ArgumentParser(
        description="Advanced Stellar Structure Plotter - Unified visualization tool for compact object analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Basic mass-radius plot with publication theme
  python stellar_plotter.py mass-radius --theme publication
  
  # Dark theme for presentations
  python stellar_plotter.py compare --theme dark --eos-types hybrid neutron_relativistic
  
  # Colorblind-friendly analysis
  python stellar_plotter.py all --theme colorblind --validate-strict
  
  # High-resolution plots for publication
  python stellar_plotter.py mass-radius --dpi 600 --output-prefix publication_
  
  # Custom styling with validation
  python stellar_plotter.py profile --file hybrid_rhoc_1.00pp16.csv --theme presentation --no-grid
  
THEMES:
  publication  - Clean, serif fonts, publication-ready (default)
  presentation - Sans-serif, larger fonts, light background
  dark         - Dark theme for presentations/screens
  colorblind   - Colorblind-friendly palette
  
VALIDATION LEVELS:
  none    - No validation (fastest)
  basic   - Basic range checks (default)
  strict  - Strict physical consistency checks
        """
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available plotting commands')
    
    # Profile command (enhanced)
    profile_parser = subparsers.add_parser('profile', help='Single stellar profile plot')
    profile_parser.add_argument('--file', required=True, help='CSV file for profile plot')
    profile_parser.add_argument('--eos-type', help='EOS type (auto-detected if not specified)')
    
    # Mass-radius command (enhanced)
    mr_parser = subparsers.add_parser('mass-radius', help='Mass-radius relation plots')
    mr_parser.add_argument('--eos-types', nargs='+', 
                          choices=['hybrid', 'neutron_relativistic', 'electron_relativistic',
                                  'electron_non_relativistic', 'neutron_non_relativistic'],
                          default=['hybrid'], help='EOS types to include')
    
    # Mass-density command (enhanced)
    md_parser = subparsers.add_parser('mass-density', help='Mass vs central density plots')
    md_parser.add_argument('--eos-types', nargs='+',
                          choices=['hybrid', 'neutron_relativistic', 'electron_relativistic',
                                  'electron_non_relativistic', 'neutron_non_relativistic'],
                          default=['hybrid'], help='EOS types to include')
    
    # Compare command (enhanced)
    compare_parser = subparsers.add_parser('compare', help='Comparative analysis plots')
    compare_parser.add_argument('--eos-types', nargs='+',
                               choices=['hybrid', 'neutron_relativistic', 'electron_relativistic',
                                       'electron_non_relativistic', 'neutron_non_relativistic'],
                               default=['hybrid', 'neutron_relativistic'], help='EOS types to compare')
    compare_parser.add_argument('--plot-types', nargs='+', choices=['mass_radius', 'mass_density'],
                               default=['mass_radius', 'mass_density'], help='Types of comparison plots')
    
    # All command (enhanced)
    all_parser = subparsers.add_parser('all', help='Generate all available plots')
    all_parser.add_argument('--eos-types', nargs='+',
                           choices=['hybrid', 'neutron_relativistic', 'electron_relativistic',
                                   'electron_non_relativistic', 'neutron_non_relativistic'],
                           default=['hybrid'], help='EOS types to include')
    
    # Global options (enhanced for Sprint 3)
    for subparser in [profile_parser, mr_parser, md_parser, compare_parser, all_parser]:
        # Theme and styling options
        subparser.add_argument('--theme', choices=['publication', 'presentation', 'dark', 'colorblind'],
                              default='publication', help='Plot theme (default: publication)')
        subparser.add_argument('--dpi', type=int, default=300, help='Output resolution (default: 300)')
        subparser.add_argument('--figure-size', nargs=2, type=float, default=[10.0, 6.0],
                      help='Figure size in inches (e.g., 3.4 2.3)')
        subparser.add_argument('--no-grid', action='store_true', help='Disable grid lines')
        subparser.add_argument('--alpha', type=float, default=0.7, help='Marker transparency (default: 0.7)')
        subparser.add_argument('--marker-size', type=float, default=30.0, help='Marker size (default: 30)')
        subparser.add_argument('--line-width', type=float, default=2.0, help='Line width (default: 2.0)')
        
        # Validation options
        subparser.add_argument('--validate', choices=['none', 'basic', 'strict'], default='basic',
                              help='Data validation level (default: basic)')
        subparser.add_argument('--validate-strict', action='store_true',
                              help='Enable strict validation (equivalent to --validate strict)')
        
        # Performance options
        subparser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing')
        subparser.add_argument('--chunk-size', type=int, default=10, help='Processing chunk size (default: 10)')
        
        # Output options
        subparser.add_argument('--config', help='YAML configuration file')
        subparser.add_argument('--output-dir', default=None, help='Output directory (default: ./plots/)')
        subparser.add_argument('--output-prefix', default='', help='Output filename prefix')
        subparser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                              default='INFO', help='Logging level (default: INFO)')
        subparser.add_argument('--stats', action='store_true', help='Display summary statistics')
        
        # Data options
        subparser.add_argument('--data-dir', default=None, help='Data directory (default: ./data/)')

        subparser.add_argument('--no-title', action='store_true')
        subparser.add_argument('--no-legend', action='store_true')
        subparser.add_argument('--connect-points', action='store_true')
        subparser.add_argument('--mark-max', action='store_true')

    return parser


def execute_plotting_command(args) -> None:
    """Execute plotting command with enhanced Sprint 3 features."""
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)
    
    try:
        # Load and customize configuration
        config_manager = ConfigManager(args.config if hasattr(args, 'config') and args.config else None)
        config = config_manager.get_config()
        
        # Apply CLI overrides to configuration
        if getattr(args, 'output_dir', None) is not None:
            config.output_directory = args.output_dir
        if getattr(args, 'data_dir', None) is not None:
            config.data_directory = args.data_dir
        if getattr(args, 'no_title', False):    config.style.show_title = False
        if getattr(args, 'no_legend', False):   config.style.show_legend = False
        if getattr(args, 'connect_points', False): config.style.connect_points = True
        if getattr(args, 'mark_max', False):    config.style.mark_max_point = True
            
        # Apply styling overrides
        if hasattr(args, 'theme'):
            config.style.theme = args.theme
        if hasattr(args, 'dpi'):
            config.style.dpi = args.dpi
        if hasattr(args, 'figure_size'):
            config.style.figure_size = tuple(args.figure_size)
        if hasattr(args, 'no_grid') and args.no_grid:
            config.style.grid = False
        if hasattr(args, 'alpha'):
            config.style.alpha = args.alpha
        if hasattr(args, 'marker_size'):
            config.style.marker_size = args.marker_size
        if hasattr(args, 'line_width'):
            config.style.line_width = args.line_width
            
        # Apply validation settings
        if hasattr(args, 'validate_strict') and args.validate_strict:
            args.validate = 'strict'
        if hasattr(args, 'validate'):
            if args.validate == 'none':
                config.validation = None
            elif args.validate == 'strict':
                # Stricter validation limits
                config.validation.min_mass_solar = 0.05
                config.validation.max_mass_solar = 3.0
                config.validation.min_radius_km = 2.0
                config.validation.max_radius_km = 50.0
                
        # Apply performance settings
        if hasattr(args, 'no_parallel') and args.no_parallel:
            config.performance.enable_parallel_loading = False
        if hasattr(args, 'chunk_size'):
            config.performance.chunk_size = args.chunk_size
        
        # Validate configuration
        if not config_manager.validate_paths():
            logger.error("Configuration validation failed")
            return
        
        # Initialize components
        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)
        
        logger.info(f"Executing command: {args.command}")
        logger.info(f"Theme: {config.style.theme}, Validation: {getattr(args, 'validate', 'basic')}")
        logging.getLogger(__name__).info(
            f"Using data dir: {config.data_directory} | output dir: {config.output_directory}"
        )

        # Execute specific command
        if args.command == 'profile':
            # Enhanced profile plotting
            if not os.path.exists(args.file):
                logger.error(f"File not found: {args.file}")
                return
            
            # Auto-detect EOS type if not provided
            eos_type = getattr(args, 'eos_type', None)
            if not eos_type:
                for eos in config.eos_patterns.keys():
                    if eos in args.file:
                        eos_type = eos
                        break
                if not eos_type:
                    logger.error("Could not auto-detect EOS type. Please specify --eos-type")
                    return
            
            model = processor.load_stellar_model(args.file, eos_type, load_profile=True)
            if not model:
                logger.error(f"Failed to load model from {args.file}")
                return
            
            output_file = f"{args.output_prefix}profile_{model.eos_type}_rhoc_{model.central_density:.2e}.png".replace('+', 'p').replace('e', 'p')
            plot_path = plotter.plot_single_profile(model, output_file)
            logger.info(f"Profile plot saved: {plot_path}")
            
        elif args.command in ['mass-radius', 'mass-density']:
            # Enhanced mass-radius and mass-density plotting
            datasets = []
            for eos_type in args.eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)
                    logger.info(f"Loaded {len(dataset.models)} models for {eos_type}")
            
            if not datasets:
                logger.error("No valid datasets loaded")
                return
            
            if args.command == 'mass-radius':
                output_file = f"{args.output_prefix}mass_radius_relation.png"
                plot_path = plotter.plot_mass_radius_relation(datasets, output_file)
            else:  # mass-density
                output_file = f"{args.output_prefix}mass_density_relation.png"
                plot_path = plotter.plot_mass_density_relation(datasets, output_file)
            
            logger.info(f"Plot saved: {plot_path}")
            
        elif args.command == 'compare':
            # Enhanced comparative analysis
            datasets = []
            for eos_type in args.eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)
            
            if not datasets:
                logger.error("No valid datasets loaded for comparison")
                return
            
            plot_types = getattr(args, 'plot_types', ['mass_radius', 'mass_density'])
            plot_paths = plotter.plot_comparative_analysis(datasets, plot_types, args.output_prefix)
            for path in plot_paths:
                logger.info(f"Comparative plot saved: {path}")
                
        elif args.command == 'all':
            # Enhanced comprehensive analysis
            datasets = []
            for eos_type in args.eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)
            
            if not datasets:
                logger.error("No valid datasets loaded")
                return
            
            # Generate all plot types
            plot_paths = []
            
            # Mass-radius relation
            mr_file = f"{args.output_prefix}comprehensive_mass_radius.png"
            plot_paths.append(plotter.plot_mass_radius_relation(datasets, mr_file))
            
            # Mass-density relation
            md_file = f"{args.output_prefix}comprehensive_mass_density.png"
            plot_paths.append(plotter.plot_mass_density_relation(datasets, md_file))
            
            # Comparative analysis
            comp_paths = plotter.plot_comparative_analysis(datasets, ['mass_radius', 'mass_density'], 
                                                         f"{args.output_prefix}comparative_")
            plot_paths.extend(comp_paths)
            
            for path in plot_paths:
                logger.info(f"Plot saved: {path}")
        
        # Display summary statistics if requested
        if hasattr(args, 'stats') and args.stats:
            if 'datasets' in locals() and datasets:
                stats = plotter.create_summary_statistics(datasets)
                logger.info("=== SUMMARY STATISTICS ===")
                for eos_type, eos_stats in stats.items():
                    logger.info(f"\n{eos_type.replace('_', ' ').title()}:")
                    logger.info(f"  Models: {eos_stats['num_models']}")
                    logger.info(f"  Mass range: {eos_stats['min_mass']:.3f} - {eos_stats['max_mass']:.3f} M☉")
                    logger.info(f"  Radius range: {eos_stats['min_radius']:.1f} - {eos_stats['max_radius']:.1f} km")
                    logger.info(f"  Density log range: {eos_stats['density_range_log']:.1f}")
        
        logger.info("Plotting completed successfully!")
        
    except Exception as e:
        logger.error(f"Error during plotting: {e}")
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