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
from scipy.interpolate import interp1d
# import yaml  # Removed for streamlined CLI
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import sys


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

_TOV_BASENAME_RE = re.compile(
    r'^tov_solution_(?P<eos>.+)_(?P<dens>(?:\d+(?:\.\d+)?(?:e[+-]?\d+)?|\d+(?:\.\d+)?pp\d+|\d+(?:\.\d+)?p-\d+))\.csv$',
    re.IGNORECASE
)


def parse_density_token(token: str) -> float:
    """Parse density token from filename into numeric value."""
    t = token.strip().lower()
    if 'pp' in t:
        base, exp = t.split('pp', 1)
        return float(base) * 10.0 ** float(exp)
    if 'p-' in t:
        base, exp = t.split('p-', 1)
        return float(base) * 10.0 ** (-float(exp))
    return float(t)


def _norm_eos_id(s: str) -> str:
    """Normalize EOS ID for consistent matching."""
    s = s.strip().lower()
    s = re.sub(r'[\s\-]+', '_', s)
    return re.sub(r'_+', '_', s).strip('_')


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
    theme: str = 'default'  # 'default', 'dark'
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
        # Map new theme names to internal configurations
        theme_mapping = {
            'default': 'publication',
            'dark': 'dark',
            # Legacy support (internal use only - silently map)
            'publication': 'publication',
            'presentation': 'publication',
            'colorblind': 'publication',
        }

        internal_theme = theme_mapping.get(theme_name, 'publication')

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
            }
        }

        return themes.get(internal_theme, themes['publication'])

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
    max_radius_km: float = 50.0       # Maximum reasonable radius (km) - limit to compact object regime
    min_density_log: float = 10.0     # Minimum log density (g/cm³)
    max_density_log: float = 20.0     # Maximum log density (g/cm³)
    skip_r_end_artifacts: bool = True       # enable by default
    r_end_km: float = 50.0                  # Changed from 1000 km to 50 km for more inclusive validation
    r_end_tolerance_km: float = 1.0         # how close to r_end_km counts as "ended at r_end"

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
            config_file: Optional path to YAML configuration file (ignored)
        """
        self.logger = logging.getLogger(__name__)
        self.config = self._load_config(config_file)

    def _load_config(self, config_file: Optional[str]) -> PlotConfig:
        """Load configuration (always returns defaults for streamlined CLI)."""
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
# EOS TABLE PARSER
# ============================================================================

DENSITY_CANDIDATES = ['log_rho','logrho','log_density','rho','density','log_n','n','number_density']
PRESSURE_CANDIDATES = ['log_p','logp','log_pressure','p','pressure']


class EOSTableParser:
    """Handles parsing of arbitrary EOS table files with robust column detection."""

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def _normalize_header(self, header: str) -> str:
        """Normalize column header for consistent matching."""
        h = header.strip().lower()
        h = h.replace(' ', '').replace('_', '')
        # Remove common unit annotations
        for unit_pattern in ['[g/cm^3]', '(g/cm^3)', '[g/cm³]', '(g/cm³)',
                           '[dyne/cm^2]', '(dyne/cm^2)', '[dyne/cm²]', '(dyne/cm²)',
                           '[g]', '(g)', '[cm^3]', '(cm^3)', '[cm³]', '(cm³)']:
            h = h.replace(unit_pattern, '')
        return h

    def _find_column(self, normalized_headers: List[str], candidates: List[str]) -> Optional[str]:
        """Find best matching column from candidates list."""
        for candidate in candidates:
            norm_candidate = self._normalize_header(candidate)
            for i, norm_header in enumerate(normalized_headers):
                if norm_candidate == norm_header or norm_candidate in norm_header:
                    return normalized_headers[i]  # Return the normalized version
        return None

    def _infer_delimiter(self, filepath: str) -> pd.DataFrame:
        """Try to infer the correct delimiter and read the file."""
        try:
            df = pd.read_csv(filepath, comment='#')
            if len(df.columns) == 1 and df.shape[0] > 0:
                cell = str(df.iloc[0, 0])
                # whitespace or tab-separated heuristic
                if bool(re.search(r'\s+', cell)) and not (',' in cell):
                    self.logger.info(f"Detected whitespace-delimited format in {filepath}")
                    df = pd.read_csv(filepath, sep=r'\s+', comment='#', engine='python')
            return df
        except Exception:
            try:
                return pd.read_csv(filepath, sep=r'\s+', comment='#', engine='python')
            except Exception as e2:
                self.logger.error(f"Failed to read {filepath}: {e2}")
                return None

    def parse_eos_file(self, filepath: str) -> Optional[Tuple[np.ndarray, np.ndarray, str]]:
        """
        Parse EOS file and return (log_rho, log_P, label).

        Args:
            filepath: Path to EOS table file

        Returns:
            Tuple of (log_density_array, log_pressure_array, label) or None if failed
        """
        df = self._infer_delimiter(filepath)
        if df is None:
            return None

        # Drop all-NaN rows
        df = df.dropna(how='all')

        if len(df) < 3:
            self.logger.error(f"Insufficient data in {filepath}: only {len(df)} valid rows")
            return None

        # Normalize headers for matching
        original_headers = list(df.columns)
        normalized_headers = [self._normalize_header(h) for h in original_headers]

        # Find density column
        density_col_norm = self._find_column(normalized_headers, DENSITY_CANDIDATES)
        if density_col_norm is None:
            self.logger.error(f"Could not find density column in {filepath}. Available columns: {original_headers}")
            return None

        # Find pressure column
        pressure_col_norm = self._find_column(normalized_headers, PRESSURE_CANDIDATES)
        if pressure_col_norm is None:
            self.logger.error(f"Could not find pressure column in {filepath}. Available columns: {original_headers}")
            return None

        # Map back to original column names
        density_col_idx = normalized_headers.index(density_col_norm)
        pressure_col_idx = normalized_headers.index(pressure_col_norm)
        density_col = original_headers[density_col_idx]
        pressure_col = original_headers[pressure_col_idx]

        # Extract data
        try:
            density_data = pd.to_numeric(df[density_col], errors='coerce').dropna()
            pressure_data = pd.to_numeric(df[pressure_col], errors='coerce').dropna()

            # Align the arrays (keep only rows where both are valid)
            valid_indices = df[density_col].notna() & df[pressure_col].notna()
            density_data = pd.to_numeric(df.loc[valid_indices, density_col])
            pressure_data = pd.to_numeric(df.loc[valid_indices, pressure_col])

            if len(density_data) < 3:
                self.logger.error(f"Insufficient valid data points in {filepath}")
                return None

        except Exception as e:
            self.logger.error(f"Failed to parse numeric data from {filepath}: {e}")
            return None

        # Auto log detection and conversion
        density_is_log = any(log_indicator in density_col.lower()
                           for log_indicator in ['log', 'lg', 'ln'])
        pressure_is_log = any(log_indicator in pressure_col.lower()
                            for log_indicator in ['log', 'lg', 'ln'])

        # Detect ln columns for conversion
        density_has_ln = density_is_log and 'ln' in density_col.lower()
        pressure_has_ln = pressure_is_log and 'ln' in pressure_col.lower()

        if density_is_log:
            log_density = density_data.values
            if density_has_ln:
                log_density = log_density / np.log(10.0)
        else:
            # Check for non-positive values
            non_positive = density_data <= 0
            if non_positive.any():
                bad_indices = non_positive[non_positive].index.tolist()
                self.logger.error(f"Non-positive density values in {filepath}: {non_positive.sum()} values, first at row {bad_indices[0]}")
                return None
            log_density = np.log10(density_data.values)

        if pressure_is_log:
            log_pressure = pressure_data.values
            if pressure_has_ln:
                log_pressure = log_pressure / np.log(10.0)
        else:
            # Check for non-positive values
            non_positive = pressure_data <= 0
            if non_positive.any():
                bad_indices = non_positive[non_positive].index.tolist()
                self.logger.error(f"Non-positive pressure values in {filepath}: {non_positive.sum()} values, first at row {bad_indices[0]}")
                return None
            log_pressure = np.log10(pressure_data.values)

        # Determine label
        label = os.path.splitext(os.path.basename(filepath))[0]

        # Check for explicit label/name columns
        for label_candidate in ['label', 'name', 'eos']:
            if label_candidate in normalized_headers:
                label_col_idx = normalized_headers.index(label_candidate)
                label_col = original_headers[label_col_idx]
                if not df[label_col].isna().all():
                    label = str(df[label_col].iloc[0])
                break

        self.logger.info(f"Parsed {filepath}: {len(log_density)} points, label='{label}'")
        return log_density, log_pressure, label


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

    def auto_discover_eos_types(self) -> Dict[str, List[str]]:
        """
        Auto-discover EOS types from tov_solution_*.csv files in data directory.
        Handles magnetic field subdirectories (b_*) by incorporating field strength into EOS type.

        Returns:
            Dictionary mapping canonical normalized EOS IDs to lists of files,
            sorted by numeric central density
        """
        # Search recursively in data directory and subdirectories
        pattern = os.path.join(self.config.data_directory, '**', 'tov_solution_*.csv')
        files = glob.glob(pattern, recursive=True)

        if not files:
            self.logger.info("No tov_solution_*.csv files found for auto-discovery")
            return {}

        discovered_eos = {}
        files_with_density = {}  # For sorting by density

        for filepath in files:
            basename = os.path.basename(filepath)
            match = _TOV_BASENAME_RE.match(basename)

            if match:
                eos_id = match.group('eos')
                dens_token = match.group('dens')

                # Check if file is in a magnetic field subdirectory (b_*)
                relative_path = os.path.relpath(filepath, self.config.data_directory)
                path_parts = relative_path.split(os.sep)

                # If file is in a subdirectory that starts with 'b_', incorporate field strength
                if len(path_parts) > 1 and path_parts[0].startswith('b_'):
                    b_value = path_parts[0]  # e.g., 'b_1e-02'
                    # Extract field value after 'b_'
                    field_str = b_value[2:]  # Remove 'b_' prefix
                    enhanced_eos_id = f"{eos_id}_b_{field_str}"
                else:
                    enhanced_eos_id = eos_id

                try:
                    density = parse_density_token(dens_token)
                    canonical_eos = _norm_eos_id(enhanced_eos_id)

                    if canonical_eos not in discovered_eos:
                        discovered_eos[canonical_eos] = []
                        files_with_density[canonical_eos] = []

                    discovered_eos[canonical_eos].append(filepath)
                    files_with_density[canonical_eos].append((density, filepath))

                except (ValueError, IndexError) as e:
                    self.logger.warning(f"Could not parse density from {basename}: {e}")
            else:
                # Only warn for files that start with tov_solution_ but fail regex
                if basename.startswith('tov_solution_'):
                    self.logger.warning(f"Could not parse tov_solution pattern from: {basename}")

        # Sort files by numeric density for each EOS type
        for canonical_eos in discovered_eos:
            files_with_density[canonical_eos].sort(key=lambda x: x[0])  # Sort by density
            discovered_eos[canonical_eos] = [filepath for _, filepath in files_with_density[canonical_eos]]
            self.logger.info(f"Auto-discovered {len(discovered_eos[canonical_eos])} files for EOS: {canonical_eos}")

        return discovered_eos

    def discover_files(self, eos_types: Optional[List[str]] = None) -> Dict[str, List[str]]:
        """
        Discover CSV files for specified EOS types using auto-discovery as primary method.

        Args:
            eos_types: List of EOS type identifiers. If None or empty, return all discovered groups.

        Returns:
            Dictionary mapping canonical EOS IDs to lists of file paths
        """
        # Primary method: Auto-discovery
        auto_discovered = self.auto_discover_eos_types()

        # If no filters specified, return all discovered groups
        if not eos_types:
            return auto_discovered

        # If primary found groups, apply fuzzy matching
        if auto_discovered:
            discovered_files = {}

            for requested_eos in eos_types:
                norm_request = _norm_eos_id(requested_eos)
                matches = []

                for canonical_eos in auto_discovered.keys():
                    # Fuzzy rule: all request tokens must be substrings of the canonical EOS
                    request_tokens = norm_request.split('_')
                    if all(token in canonical_eos for token in request_tokens if token):
                        matches.append(canonical_eos)

                if matches:
                    for matched_eos in matches:
                        discovered_files[matched_eos] = auto_discovered[matched_eos]
                        self.logger.info(f"Fuzzy match: '{requested_eos}' -> {matched_eos} ({len(auto_discovered[matched_eos])} files)")
                else:
                    self.logger.warning(f"No auto-discovered EOS matches '{requested_eos}'. Available: {list(auto_discovered.keys())}")

            return discovered_files

        # Fallback method: Use hardcoded patterns if primary scan found nothing
        self.logger.info("No auto-discovery results, falling back to hardcoded patterns")
        fallback_files = {}

        for eos_type in eos_types:
            if eos_type in self.config.eos_patterns:
                # Try recursive search for hardcoded patterns
                pattern = os.path.join(self.config.data_directory, '**',
                                     os.path.basename(self.config.eos_patterns[eos_type]))
                files = glob.glob(pattern, recursive=True)

                if files:
                    canonical_key = _norm_eos_id(eos_type)
                    fallback_files[canonical_key] = sorted(files)
                    self.logger.info(f"Hardcoded pattern match: {eos_type} -> {canonical_key} ({len(files)} files)")
                else:
                    self.logger.warning(f"No files found for {eos_type} with pattern: {pattern}")
            else:
                self.logger.warning(f"Unknown EOS type: {eos_type}")

        return fallback_files

    def parse_central_density(self, filename: str, eos_type: str) -> Optional[float]:
        """
        Parse central density from filename.

        Args:
            filename: Path to CSV file
            eos_type: EOS type identifier (for compatibility checking)

        Returns:
            Central density in g/cm³, or None if parsing fails
        """
        base = os.path.basename(filename)

        try:
            # Primary: Try standard tov_solution pattern
            match = _TOV_BASENAME_RE.match(base)
            if match:
                dens_token = match.group('dens')
                return parse_density_token(dens_token)

            # Fallback: Legacy hardcoded patterns for backwards compatibility
            if eos_type == 'hybrid':
                # Pattern: tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv
                legacy_match = re.search(r'polytrope_(.*)\.csv', base)
                if legacy_match:
                    rho_str = legacy_match.group(1)  # e.g. '1.00e+16'
                    return parse_density_token(rho_str)

            else:
                # Pattern: neutron_relativistic_rhoc_5.00pp18.csv
                legacy_match = re.search(rf'{eos_type}_rhoc_(.*)\.csv', base)
                if legacy_match:
                    rho_str = legacy_match.group(1)  # e.g. '5.00pp18'
                    return parse_density_token(rho_str)

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
        groups = self.discover_files([eos_type])
        files = []
        if groups:
            # collect all matched groups (fuzzy may return >1)
            for _canonical, _files in groups.items():
                files.extend(_files)
        files = sorted(set(files))
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

        # Colorblind theme no longer supported via CLI (legacy code removed)

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

    def plot_eos_table(self, eos_files: List[str], output_file: str = "eos_comparison.png") -> str:
        """
        Create EOS comparison plot from arbitrary EOS table files.

        Args:
            eos_files: List of paths to EOS table files
            output_file: Output filename

        Returns:
            Path to saved plot file
        """
        parser = EOSTableParser()
        eos_data = []

        for filepath in eos_files:
            result = parser.parse_eos_file(filepath)
            if result is not None:
                log_density, log_pressure, label = result
                eos_data.append((log_density, log_pressure, label))
            else:
                self.logger.warning(f"Failed to parse {filepath}")

        if not eos_data:
            self.logger.error("No valid EOS files could be parsed")
            return None

        # Create the plot
        fig, ax = plt.subplots(figsize=self.config.style.figure_size)

        for log_density, log_pressure, label in eos_data:
            ax.scatter(log_density, log_pressure,
                      s=10,  # Small marker size for dense data
                      marker='o',
                      alpha=0.6,
                      label=label)

        # Labels and formatting
        ax.set_xlabel(r'$\log_{10}\rho$ [g cm$^{-3}$]')
        ax.set_ylabel(r'$\log_{10}P$ [dyne cm$^{-2}$]')

        # Apply current theme styling
        ax.tick_params(direction='in', which='both', top=True, right=True, length=4)
        ax.tick_params(which='minor', length=2)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)

        if self.config.style.show_title:
            ax.set_title('Equation of State', pad=10)

        if self.config.style.show_legend and len(eos_data) > 1:
            ax.legend(loc=self.config.style.legend_location)

        plt.tight_layout()

        # Save both PNG and PDF
        output_path = os.path.join(self.config.output_directory, output_file)
        fig.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')

        pdf_path = os.path.splitext(output_path)[0] + ".pdf"
        fig.savefig(pdf_path, bbox_inches='tight')

        plt.close(fig)

        self.logger.info(f"Saved EOS comparison plot: {output_path}")
        self.logger.info(f"Saved EOS comparison plot: {pdf_path}")

        return output_path

    # ========================================================================
    # MAGNETIC FIELD STRATIFIED PLOTTING
    # ========================================================================

    def _list_magnetic_field_dirs(self, root_dir: str) -> List[Tuple[float, str]]:
        """
        Scan for b_* subdirectories and extract magnetic field values.

        Args:
            root_dir: Root directory to search for b_* subdirectories

        Returns:
            List of (b_value, directory_path) tuples, sorted by b_value
        """
        pairs = []
        pattern = os.path.join(root_dir, "b_*")

        for d in sorted(glob.glob(pattern)):
            base = os.path.basename(d)
            match = re.match(r"b_([0-9.eE+-]+)$", base)
            if not match:
                continue
            try:
                b = float(match.group(1))
                pairs.append((b, d))
            except ValueError:
                self.logger.warning(f"Could not parse b value from directory: {base}")
                continue

        return sorted(pairs, key=lambda x: x[0])

    def _read_surface_values(self, csv_path: str, r_end_cm: float) -> Optional[Tuple[float, float]]:
        """
        Extract surface mass and radius from a TOV solution CSV file.

        Args:
            csv_path: Path to CSV file
            r_end_cm: Computational boundary radius in cm (for artifact detection)

        Returns:
            Tuple of (M_solar, R_km) or None if invalid
        """
        SOLAR_MASS_G = 1.989e33
        CM_TO_KM = 1e-5

        try:
            df = pd.read_csv(csv_path)

            # Check required columns
            needed = {'log_m[g]', 'log_r[cm]'}
            if not needed.issubset(df.columns):
                return None

            # Optional pressure monotonicity check
            if 'log_P[dyne/cm^2]' in df.columns and len(df) > 3:
                lp = df['log_P[dyne/cm^2]'].to_numpy()
                # Pressure should decrease monotonically
                if not np.all(np.diff(lp) <= 1e-6):
                    return None
                # Last pressure should be close to minimum
                if not np.isclose(lp[-1], lp.min(), atol=5e-4):
                    return None

            # Extract surface values
            m_sol = 10.0**df['log_m[g]'].iloc[-1] / SOLAR_MASS_G
            r_km = 10.0**df['log_r[cm]'].iloc[-1] * CM_TO_KM

            # Skip models that hit computational boundary
            if r_end_cm is not None:
                r_end_km = r_end_cm * CM_TO_KM
                if (r_km > 0.98*r_end_km) and np.isclose(r_km, r_end_km, rtol=1e-3, atol=0.5):
                    return None

            return m_sol, r_km

        except Exception as e:
            self.logger.debug(f"Failed to read {csv_path}: {e}")
            return None

    def _collect_mr_family(self, b_dir: str, r_end_cm: float, file_list: Optional[List[str]] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Collect all valid M-R pairs from a magnetic field directory.

        Args:
            b_dir: Directory containing CSV files for a specific b value
            r_end_cm: Computational boundary radius in cm
            file_list: Optional list of specific files to use (default: None = use all CSV files)

        Returns:
            Tuple of (R_km_sorted, M_solar_sorted) arrays
        """
        all_pairs = []

        # Use provided file list if available, otherwise glob all CSV files
        if file_list is not None:
            csv_files = file_list
        else:
            csv_files = glob.glob(os.path.join(b_dir, "*.csv"))

        for csv in csv_files:
            mr = self._read_surface_values(csv, r_end_cm=r_end_cm)
            if mr is not None:
                all_pairs.append(mr)

        if not all_pairs:
            return np.array([]), np.array([])

        M, R = zip(*all_pairs)
        M = np.array(M)
        R = np.array(R)

        # Sort by radius for smooth plotting
        idx = np.argsort(R)
        return R[idx], M[idx]

    def plot_mass_radius_by_field(
        self,
        data_dir: str,
        r_end_cm: float = 1.0e8,
        output_file: str = "mr_by_b.png",
        eos_types: Optional[List[str]] = None
    ) -> str:
        """
        Create mass-radius plot stratified by magnetic field strength.

        Generates a two-panel figure:
        - Top panel: M-R curves for different b values
        - Bottom panel: Residuals (ΔM) relative to reference field

        Args:
            data_dir: Root directory containing b_* subdirectories
            r_end_cm: Computational boundary radius in cm (default: 1e8 cm = 1000 km)
            output_file: Output filename (default: mr_by_b.png)
            eos_types: Optional list of EOS types to filter (default: None = all files)

        Returns:
            Path to saved plot file
        """
        # Discover b_* directories
        b_dirs = self._list_magnetic_field_dirs(data_dir)
        if not b_dirs:
            self.logger.error(f"No b_* directories found in {data_dir}")
            return None

        # Apply EOS type filtering if requested
        filtered_files_by_b = {}
        if eos_types is not None and len(eos_types) > 0:
            # Use EOSDataProcessor for fuzzy matching across all b_* directories
            processor = EOSDataProcessor(self.config)

            # Discover files recursively and apply fuzzy matching
            discovered_groups = processor.discover_files(eos_types)

            # Group discovered files by their b_* directory
            for canonical_eos, file_list in discovered_groups.items():
                for filepath in file_list:
                    # Check which b_* directory this file belongs to
                    for b, b_path in b_dirs:
                        if filepath.startswith(b_path):
                            if b not in filtered_files_by_b:
                                filtered_files_by_b[b] = []
                            filtered_files_by_b[b].append(filepath)
                            break

            if not filtered_files_by_b:
                self.logger.warning(f"No files matching EOS types {eos_types} found in any b_* directory")
                return None

            self.logger.info(f"EOS filtering: selected {sum(len(files) for files in filtered_files_by_b.values())} files across {len(filtered_files_by_b)} field strengths")

        # Determine radius cap for filtering (matching plot_mass_radius_relation behavior)
        if getattr(self.config, "validation", None) and hasattr(self.config.validation, "max_radius_km"):
            radius_cap = float(self.config.validation.max_radius_km)
            self.logger.info(f"Applying radius filter: max_radius = {radius_cap} km")
        else:
            radius_cap = float("inf")

        def apply_radius_filter(radii: np.ndarray, masses: np.ndarray, context: str, log_drops: bool = True) -> Tuple[np.ndarray, np.ndarray, int]:
            """Apply radius cap to radius/mass points, optionally logging dropped counts."""
            if not math.isfinite(radius_cap) or radii.size == 0:
                return radii, masses, 0
            mask = radii <= radius_cap
            kept = int(mask.sum())
            dropped = int(mask.size - kept)
            if dropped > 0 and log_drops:
                self.logger.info(f"{context}: dropped {dropped} points with R > {radius_cap} km")
            return radii[mask], masses[mask], dropped

        # Create two-panel figure with gridspec
        fig = plt.figure(figsize=(7.5, 7.0), constrained_layout=True)
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
        ax = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax)

        ref_b = None
        ref_curve = None
        summary = []

        # Plot M-R curves for each b value
        for b, path in b_dirs:
            # Use filtered file list if available, otherwise use all files
            if filtered_files_by_b:
                if b in filtered_files_by_b:
                    file_list = sorted(filtered_files_by_b[b])
                else:
                    # Skip this b value if no matching files after EOS filtering
                    continue
            else:
                file_list = None  # Will use all CSV files in the directory

            files = sorted(glob.glob(os.path.join(path, "*.csv"))) if file_list is None else file_list
            R_km, M_solar = self._collect_mr_family(path, r_end_cm=r_end_cm, file_list=file_list)
            R_km, M_solar, _ = apply_radius_filter(R_km, M_solar, path)
            used = len(R_km)

            summary.append((
                b, len(files), used,
                (float(np.min(R_km)) if used else np.nan,
                 float(np.max(R_km)) if used else np.nan)
            ))

            if used == 0:
                self.logger.warning(f"{path}: {len(files)} files, 0 usable (likely r_end/surface filter)")
                continue

            label = f"b = {b:g}"
            ax.scatter(R_km, M_solar, s=9.0, alpha=0.9, zorder=2, label=label)

            if self.config.style.grid:
                ax.grid(True, alpha=self.config.style.grid_alpha, zorder=0)

            # Pick reference as smallest b
            if ref_curve is None or (ref_b is not None and b < ref_b):
                ref_b, ref_curve = b, (R_km, M_solar)

        # Compute residuals vs reference
        if ref_curve is not None:
            Rref, Mref = ref_curve
            fMref = interp1d(Rref, Mref, bounds_error=False, fill_value=np.nan)

            for b, path in b_dirs:
                # Use filtered file list if available
                if filtered_files_by_b:
                    if b not in filtered_files_by_b:
                        continue
                    file_list = sorted(filtered_files_by_b[b])
                else:
                    file_list = None

                R_km, M_solar = self._collect_mr_family(path, r_end_cm=r_end_cm, file_list=file_list)
                R_km, M_solar, _ = apply_radius_filter(R_km, M_solar, path, log_drops=False)
                if R_km.size == 0 or b == ref_b:
                    continue

                # Compute ΔM at overlapping radii
                Rmin = max(np.min(R_km), np.min(Rref))
                Rmax = min(np.max(R_km), np.max(Rref))
                mask = (R_km >= Rmin) & (R_km <= Rmax)
                dM = M_solar[mask] - fMref(R_km[mask])
                ax2.plot(R_km[mask], dM, lw=1.2, label=f"ΔM (b={b:g} − {ref_b:g})")

        # Format top panel
        if self.config.style.show_title:
            ax.set_title(r"Mass–Radius curves colored by $b = B/B_{\rm Q}$")
        ax.set_xlabel("Radius R (km)")
        ax.set_ylabel(r"Gravitational Mass M ($M_\odot$)")
        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)
        if self.config.style.show_legend:
            ax.legend(frameon=True, fontsize=9)

        # Format bottom panel
        ax2.axhline(0, ls=":", lw=1)
        ax2.set_xlabel("Radius R (km)")
        ax2.set_ylabel(r"ΔM ($M_\odot$)")
        if self.config.style.grid:
            ax2.grid(True, alpha=self.config.style.grid_alpha)
        if math.isfinite(radius_cap):
            ax.set_xlim(0.0, radius_cap)
            ax2.set_xlim(0.0, radius_cap)

        # Save plot
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi)
        plt.close(fig)

        self.logger.info(f"Saved magnetic field stratified M-R plot: {output_path}")

        # Log summary statistics
        self.logger.info("\nPer-b summary (files, used, R-range [km]):")
        for b, nfiles, used, Rrng in summary:
            rmin, rmax = Rrng
            self.logger.info(
                f"  b={b:g}: files={nfiles:3d}, used={used:3d}, "
                f"R∈[{'' if np.isnan(rmin) else f'{rmin:.1f}'}, "
                f"{'' if np.isnan(rmax) else f'{rmax:.1f}'}]"
            )

        return output_path

    def plot_delta_radius_vs_b(
        self,
        data_dir: str,
        target_masses: List[float],
        ref_b_tag: Optional[str] = None,
        r_end_cm: float = 1.0e8,
        output_file: str = "delta_radius_vs_b.png",
        eos_types: Optional[List[str]] = None
    ) -> str:
        """
        Compute and plot ΔR(b; M*) = R_b(M*) - R_ref(M*) at fixed masses.

        Args:
            data_dir: Root directory containing b_* subdirectories
            target_masses: List of target masses in M_sun
            ref_b_tag: Optional explicit reference family tag (e.g., 'b_0e00')
            r_end_cm: Computational boundary radius in cm
            output_file: Output filename
            eos_types: Optional list of EOS types to filter

        Returns:
            Path to saved plot file
        """
        # Discover b_* directories
        b_dirs = self._list_magnetic_field_dirs(data_dir)
        if not b_dirs:
            self.logger.error(f"No b_* directories found in {data_dir}")
            return None

        # Apply EOS type filtering if requested (reuse logic from plot_mass_radius_by_field)
        filtered_files_by_b = {}
        if eos_types is not None and len(eos_types) > 0:
            processor = EOSDataProcessor(self.config)
            discovered_groups = processor.discover_files(eos_types)

            for canonical_eos, file_list in discovered_groups.items():
                for filepath in file_list:
                    for b, b_path in b_dirs:
                        if filepath.startswith(b_path):
                            if b not in filtered_files_by_b:
                                filtered_files_by_b[b] = []
                            filtered_files_by_b[b].append(filepath)
                            break

            if not filtered_files_by_b:
                self.logger.warning(f"No files matching EOS types {eos_types} found in any b_* directory")
                return None

        # Apply radius filtering
        if getattr(self.config, "validation", None) and hasattr(self.config.validation, "max_radius_km"):
            radius_cap = float(self.config.validation.max_radius_km)
        else:
            radius_cap = float("inf")

        # Collect MR data for each family
        families = {}  # b_value -> {'M': array, 'R': array, 'tag': str, 'R_interp': callable}

        for b, path in b_dirs:
            # Use filtered file list if available
            if filtered_files_by_b:
                if b not in filtered_files_by_b:
                    continue
                file_list = sorted(filtered_files_by_b[b])
            else:
                file_list = None

            # Collect MR points (reuse existing helper)
            R_km, M_solar = self._collect_mr_family(path, r_end_cm=r_end_cm, file_list=file_list)

            # Apply radius filter
            if R_km.size > 0 and math.isfinite(radius_cap):
                mask = R_km <= radius_cap
                R_km, M_solar = R_km[mask], M_solar[mask]

            if R_km.size < 2:
                self.logger.warning(f"[warn] family {os.path.basename(path)} unusable for ΔR (insufficient stable points)")
                continue

            # Sort by mass and extract stable branch (monotonically increasing M)
            idx = np.argsort(M_solar)
            M_sorted, R_sorted = M_solar[idx], R_km[idx]

            # Keep only strictly increasing mass (stable branch up to M_max)
            keep_mask = np.concatenate([[True], np.diff(M_sorted) > 0])
            M_stable = M_sorted[keep_mask]
            R_stable = R_sorted[keep_mask]

            if M_stable.size < 2:
                self.logger.warning(f"[warn] family {os.path.basename(path)} unusable for ΔR (no stable monotone segment)")
                continue

            # Build R(M) interpolator
            R_interp = interp1d(M_stable, R_stable, kind='linear', bounds_error=False, fill_value=np.nan)

            tag = os.path.basename(path)
            families[b] = {
                'M': M_stable,
                'R': R_stable,
                'tag': tag,
                'R_interp': R_interp,
                'M_min': float(np.min(M_stable)),
                'M_max': float(np.max(M_stable))
            }

        if len(families) < 2:
            self.logger.error("Need at least 2 usable families for ΔR analysis")
            return None

        # Pick reference
        ref_b = None
        if ref_b_tag is not None:
            # User specified explicit reference
            for b, fam in families.items():
                if fam['tag'] == ref_b_tag:
                    ref_b = b
                    break
            if ref_b is None:
                self.logger.error(f"Specified reference {ref_b_tag} not found in usable families")
                return None
        else:
            # Auto-pick: b=0 if exists, else smallest b
            b_zero_candidates = [b for b in families.keys() if abs(b) < 1e-12]
            if b_zero_candidates:
                ref_b = b_zero_candidates[0]
            else:
                ref_b = min(families.keys())
                self.logger.warning(f"WARNING: using smallest b={ref_b:g} as reference (no b=0 present)")

        ref_interp = families[ref_b]['R_interp']

        # Compute mass overlap
        M_min_overlap = max(fam['M_min'] for fam in families.values())
        M_max_overlap = min(fam['M_max'] for fam in families.values())

        if M_min_overlap >= M_max_overlap:
            self.logger.error(f"No mass overlap between families: [{M_min_overlap:.2f}, {M_max_overlap:.2f}]")
            return None

        # Compute ΔR for each target mass
        results = {}  # M_target -> list of (b, ΔR) tuples
        summary = []

        for M_target in target_masses:
            if M_target < M_min_overlap or M_target > M_max_overlap:
                self.logger.warning(f"Target mass {M_target:.2f} Msun not in overlap [{M_min_overlap:.2f}, {M_max_overlap:.2f}]; skipping")
                continue

            R_ref = ref_interp(M_target)
            if np.isnan(R_ref):
                self.logger.warning(f"Reference interpolation failed for M={M_target:.2f}; skipping")
                continue

            delta_R_list = []
            for b, fam in families.items():
                R_b = fam['R_interp'](M_target)
                if np.isnan(R_b):
                    continue
                delta_R = R_b - R_ref
                delta_R_list.append((b, delta_R))

            if len(delta_R_list) < 2:
                self.logger.warning(f"Insufficient data points for M={M_target:.2f}; skipping")
                continue

            results[M_target] = sorted(delta_R_list, key=lambda x: x[0])

            # Summary statistics
            delta_R_values = [dr for _, dr in delta_R_list]
            summary.append({
                'M': M_target,
                'overlap': (M_min_overlap, M_max_overlap),
                'N_families': len(delta_R_list),
                'dR_min': min(delta_R_values),
                'dR_max': max(delta_R_values)
            })

        if not results:
            self.logger.error("No valid results to plot")
            return None

        # Create plot
        fig, ax = plt.subplots(figsize=self.config.style.figure_size)

        for M_target in sorted(results.keys()):
            b_values = [b for b, _ in results[M_target]]
            dR_values = [dr for _, dr in results[M_target]]
            ax.scatter(b_values, dR_values, s=50, label=f'M = {M_target:.2f} M$_\\odot$', alpha=0.8)

        # Use symlog to handle b=0 (log scale can't plot zero)
        # Linear in range [-linthresh, +linthresh], log outside
        b_all = [b for result_list in results.values() for b, _ in result_list]
        b_nonzero = [b for b in b_all if abs(b) > 1e-12]
        if b_nonzero:
            linthresh = min(b_nonzero) * 0.1  # Linear threshold = 10% of smallest non-zero b
        else:
            linthresh = 1e-5  # Fallback

        ax.set_xscale('symlog', linthresh=linthresh)
        ax.axhline(0, ls=':', lw=1, color='gray', alpha=0.7)
        ax.set_xlabel(r'$b = B/B_{\rm Q}$')
        ax.set_ylabel(r'$\Delta R$ (km) relative to reference')

        if self.config.style.show_title:
            ax.set_title('Radius shift at fixed mass')

        if self.config.style.grid:
            ax.grid(True, alpha=self.config.style.grid_alpha)

        if self.config.style.show_legend:
            ax.legend(loc=self.config.style.legend_location)

        plt.tight_layout()

        # Save plot
        output_path = os.path.join(self.config.output_directory, output_file)
        plt.savefig(output_path, dpi=self.config.style.dpi, bbox_inches='tight')
        plt.close(fig)

        self.logger.info(f"Saved ΔR vs b plot: {output_path}")

        # Print summary table
        print("\n" + "="*80)
        print(f"{'Target M*':<12} {'Overlap Range':<20} {'N_families':<12} {'ΔR range [km]':<20}")
        print("="*80)
        for s in summary:
            overlap_str = f"[{s['overlap'][0]:.2f}, {s['overlap'][1]:.2f}]"
            dR_range_str = f"[{s['dR_min']:+.3f}, {s['dR_max']:+.3f}]"
            print(f"{s['M']:<12.2f} {overlap_str:<20} {s['N_families']:<12} {dR_range_str:<20}")
        print("="*80 + "\n")

        return output_path


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def create_argument_parser() -> argparse.ArgumentParser:
    """Create streamlined command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Stellar Structure Plotter - Visualization tool for compact object analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Plot all discovered EOS types
  python stellar_plotter.py mr
  python stellar_plotter.py md
  python stellar_plotter.py both

  # Plot specific EOS types
  python stellar_plotter.py mr --eos-types hybrid magnetic_bps

  # Limit radius range for neutron stars
  python stellar_plotter.py mr --max-radius 20

  # Include white dwarfs (larger radii)
  python stellar_plotter.py mr --max-radius 100

  # Dark theme for presentations
  python stellar_plotter.py mr --theme dark

  # High-resolution plots with custom radius limit
  python stellar_plotter.py both --dpi 600 --max-radius 30

THEMES:
  default - Clean, publication-ready (default)
  dark    - Dark theme for presentations/screens

RADIUS FILTERING:
  --max-radius controls both data filtering and plot scaling
  Default: 50 km (compact object regime)
  Typical values: 15 km (neutron stars), 100 km (white dwarfs)
        """
    )

    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available plotting commands')

    # Core commands with aliases
    # Mass-radius (mr alias)
    mr_parser = subparsers.add_parser('mr', help='Mass-radius relation plots')
    mr_parser.add_argument('--eos-types', nargs='*', default=None,
                          help='EOS types to include (default: auto-discover all)')

    # Mass-density (md alias)
    md_parser = subparsers.add_parser('md', help='Mass vs central density plots')
    md_parser.add_argument('--eos-types', nargs='*', default=None,
                          help='EOS types to include (default: auto-discover all)')

    # Both plots
    both_parser = subparsers.add_parser('both', help='Generate both mass-radius and mass-density plots')
    both_parser.add_argument('--eos-types', nargs='*', default=None,
                            help='EOS types to include (default: auto-discover all)')

    # Magnetic field stratified M-R plots
    mrb_parser = subparsers.add_parser('mr-by-b', help='Mass-radius plots stratified by magnetic field strength')
    mrb_parser.add_argument('--eos-types', nargs='*', default=None,
                           help='EOS types to include (default: auto-discover all)')
    mrb_parser.add_argument('--r-end-cm', type=float, default=1.0e8,
                           help='Boundary cutoff in cm (default: 1e8 cm = 1000 km)')
    mrb_parser.add_argument('--output', default='mr_by_b.png',
                           help='Output filename (default: mr_by_b.png)')

    # Delta-R vs b plot
    deltaR_parser = subparsers.add_parser('deltaR', help='Radius shift at fixed mass vs magnetic field')
    deltaR_parser.add_argument('--mr-target-masses', type=str, default='1.4',
                              help='Comma/space-separated target masses in M_sun (default: 1.4)')
    deltaR_parser.add_argument('--ref-b', type=str, default=None,
                              help='Explicit reference family tag (e.g., b_0e00). Auto-picks b=0 if not specified.')
    deltaR_parser.add_argument('--eos-types', nargs='*', default=None,
                              help='EOS types to include (default: auto-discover all)')
    deltaR_parser.add_argument('--r-end-cm', type=float, default=1.0e8,
                              help='Boundary cutoff in cm (default: 1e8 cm = 1000 km)')
    deltaR_parser.add_argument('--output', type=str, default='delta_radius_vs_b.png',
                              help='Output filename (default: delta_radius_vs_b.png)')

    # Legacy aliases (hidden from help)
    mr_legacy = subparsers.add_parser('mass-radius', help=argparse.SUPPRESS)
    mr_legacy.add_argument('--eos-types', nargs='*', default=None)

    md_legacy = subparsers.add_parser('mass-density', help=argparse.SUPPRESS)
    md_legacy.add_argument('--eos-types', nargs='*', default=None)

    # Profile command
    profile_parser = subparsers.add_parser('profile', help='Single stellar profile plot')
    profile_parser.add_argument('--file', required=True, help='CSV file for profile plot')
    profile_parser.add_argument('--eos-type', help='EOS type (auto-detected if not specified)')

    # EOS command
    eos_parser = subparsers.add_parser('eos', help='Plot equation of state tables')
    eos_parser.add_argument('--files', nargs='+', required=True, help='EOS table files to plot')
    eos_parser.add_argument('--output', default='eos_comparison.png', help='Output filename (default: eos_comparison.png)')

    # List command
    list_parser = subparsers.add_parser('list', help='List auto-discovery inventory')

    # Essential options only
    for subparser in [mr_parser, md_parser, both_parser, mrb_parser, deltaR_parser, mr_legacy, md_legacy, profile_parser, eos_parser, list_parser]:
        subparser.add_argument('--theme', choices=['default', 'dark'], default='default',
                              help='Plot theme (default: default)')
        subparser.add_argument('--dpi', type=int, default=300, help='Output resolution (default: 300)')
        subparser.add_argument('--max-radius', type=float, default=None,
                              help='Maximum radius for filtering/plotting in km (default: 50.0)')
        subparser.add_argument('--output-dir', default=None, help='Output directory (default: ./plots/)')
        subparser.add_argument('--data-dir', default=None, help='Data directory (default: ./data/)')
        subparser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                              default='INFO', help='Logging level (default: INFO)')

    return parser


def data_inventory(processor: 'EOSDataProcessor', root: str) -> None:
    """Print formatted inventory of discovered TOV solutions."""
    logger = logging.getLogger(__name__)
    auto_discovered = processor.auto_discover_eos_types()
    if auto_discovered:
        logger.info(f"TOV Solutions Found in {root}:")
        logger.info("━" * 42)
        for eos_id in sorted(auto_discovered.keys()):
            files = auto_discovered[eos_id]

            # Get density range by parsing filenames
            densities = []
            for filepath in files:
                density = processor.parse_central_density(filepath, eos_id)
                if density is not None:
                    densities.append(density)

            if densities:
                min_density = min(densities)
                max_density = max(densities)
                logger.info(f"  {eos_id}")
                logger.info(f"    Files: {len(files)}")
                logger.info(f"    ρc range: {min_density:.2e} → {max_density:.2e} g/cm³")
            else:
                logger.info(f"  {eos_id}")
                logger.info(f"    Files: {len(files)} (density parsing failed)")
    else:
        logger.info(f"No TOV solutions found in {root}")


def execute_plotting_command(args) -> None:
    """Execute streamlined plotting commands."""
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)

    try:
        # Load and customize configuration
        config_manager = ConfigManager()
        config = config_manager.get_config()

        # Apply CLI overrides (essential options only)
        if args.output_dir is not None:
            config.output_directory = args.output_dir
        if args.data_dir is not None:
            config.data_directory = args.data_dir
        if hasattr(args, 'theme'):
            config.style.theme = args.theme
        if hasattr(args, 'dpi'):
            config.style.dpi = args.dpi
        if hasattr(args, 'max_radius') and args.max_radius is not None:
            config.validation.max_radius_km = args.max_radius

        # Validate configuration
        if not config_manager.validate_paths():
            logger.error("Configuration validation failed")
            return

        # Initialize components
        processor = EOSDataProcessor(config)
        plotter = StellarPlotter(config)

        # Helper function to get EOS list with auto-discovery
        def get_eos_list(eos_types_arg):
            if eos_types_arg is None or len(eos_types_arg) == 0:
                # Auto-discover all EOS
                groups = processor.discover_files(None)
                return sorted(groups.keys())
            else:
                return eos_types_arg

        # Execute specific command
        if args.command == 'profile':
            # Profile plotting
            if not os.path.exists(args.file):
                logger.error(f"File not found: {args.file}")
                return

            # Auto-detect EOS type if not provided
            eos_type = getattr(args, 'eos_type', None)
            if not eos_type:
                # Try to auto-detect from filename
                for potential_eos in processor.auto_discover_eos_types().keys():
                    if potential_eos in args.file:
                        eos_type = potential_eos
                        break
                if not eos_type:
                    logger.error("Could not auto-detect EOS type. Please specify --eos-type")
                    return

            model = processor.load_stellar_model(args.file, eos_type, load_profile=True)
            if not model:
                logger.error(f"Failed to load model from {args.file}")
                return

            output_file = f"profile_{model.eos_type}_rhoc_{model.central_density:.2e}.png".replace('+', 'p').replace('e', 'p')
            plot_path = plotter.plot_single_profile(model, output_file)
            logger.info(f"✓ Saved: {plot_path}")

        elif args.command in ['mr', 'mass-radius']:
            # Mass-radius plotting with auto-discovery
            eos_types = get_eos_list(args.eos_types)
            datasets = []
            for eos_type in eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)

            if not datasets:
                logger.error("No valid datasets loaded")
                return

            plot_path = plotter.plot_mass_radius_relation(datasets, "mass_radius.png")
            logger.info(f"✓ Saved: {plot_path}")

        elif args.command in ['md', 'mass-density']:
            # Mass-density plotting with auto-discovery
            eos_types = get_eos_list(args.eos_types)
            datasets = []
            for eos_type in eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)

            if not datasets:
                logger.error("No valid datasets loaded")
                return

            plot_path = plotter.plot_mass_density_relation(datasets, "mass_density.png")
            logger.info(f"✓ Saved: {plot_path}")

        elif args.command == 'both':
            # Both mass-radius and mass-density plots
            eos_types = get_eos_list(args.eos_types)
            datasets = []
            for eos_type in eos_types:
                dataset = processor.load_eos_dataset(eos_type)
                if dataset:
                    datasets.append(dataset)

            if not datasets:
                logger.error("No valid datasets loaded")
                return

            # Generate both plots
            mr_path = plotter.plot_mass_radius_relation(datasets, "mass_radius.png")
            md_path = plotter.plot_mass_density_relation(datasets, "mass_density.png")
            logger.info(f"✓ Saved: {mr_path}")
            logger.info(f"✓ Saved: {md_path}")

        elif args.command == 'eos':
            # EOS table plotting
            plot_path = plotter.plot_eos_table(args.files, args.output)
            if plot_path:
                logger.info(f"✓ Saved: {plot_path}")
            else:
                logger.error("Failed to create EOS comparison plot")

        elif args.command == 'mr-by-b':
            # Magnetic field stratified M-R plotting
            data_dir = args.data_dir if args.data_dir else config.data_directory
            eos_types = args.eos_types if hasattr(args, 'eos_types') else None
            plot_path = plotter.plot_mass_radius_by_field(
                data_dir=data_dir,
                r_end_cm=args.r_end_cm,
                output_file=args.output,
                eos_types=eos_types
            )
            if plot_path:
                logger.info(f"✓ Saved: {plot_path}")
            else:
                logger.error("Failed to create magnetic field stratified M-R plot")

        elif args.command == 'deltaR':
            # Delta-R vs b plotting
            # Parse target masses from comma/space-separated string
            mass_str = args.mr_target_masses
            target_masses = [float(x) for x in re.split(r'[,\s]+', mass_str.strip()) if x]

            if not target_masses:
                logger.error("No valid target masses provided")
                return

            data_dir = args.data_dir if args.data_dir else config.data_directory
            eos_types = args.eos_types if hasattr(args, 'eos_types') else None
            ref_b_tag = args.ref_b if hasattr(args, 'ref_b') else None
            r_end_cm = args.r_end_cm if hasattr(args, 'r_end_cm') else 1.0e8
            output_file = args.output if hasattr(args, 'output') else 'delta_radius_vs_b.png'

            plot_path = plotter.plot_delta_radius_vs_b(
                data_dir=data_dir,
                target_masses=target_masses,
                ref_b_tag=ref_b_tag,
                r_end_cm=r_end_cm,
                output_file=output_file,
                eos_types=eos_types
            )
            if plot_path:
                logger.info(f"✓ Saved: {plot_path}")
            else:
                logger.error("Failed to create ΔR vs b plot")

        elif args.command == 'list':
            # Auto-discovery inventory using helper
            root = os.path.abspath(config.data_directory)
            data_inventory(processor, root)

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
        print("No command specified. Available commands:")
        print("  mr | md | both | mr-by-b | deltaR | eos | profile | list")
        print("Try --help for usage details.")


if __name__ == "__main__":
    main()
