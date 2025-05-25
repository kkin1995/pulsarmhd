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
# MAIN ENTRY POINT (PLACEHOLDER)
# ============================================================================

def main():
    """Main entry point - to be expanded in Sprint 2."""
    setup_logging()
    logger = logging.getLogger(__name__)
    
    logger.info("Stellar Plotter v1.0 - Sprint 1 Foundation")
    logger.info("Core architecture initialized successfully")
    
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
    else:
        logger.error("Configuration validation failed")


if __name__ == "__main__":
    main() 