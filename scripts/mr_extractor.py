#!/usr/bin/env python3
"""
TOV Mass-Radius Data Extractor

Extracts total mass M and total radius R from TOV solution CSV files
and saves them to a structured data file for analysis.

Based on the stellar plotter extraction logic.
Author: Extracted from stellar_plotter.py
"""

import pandas as pd
import numpy as np
import glob
import os
import re
import math
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class TOVMassRadiusExtractor:
    """Extract mass and radius data from TOV solution CSV files."""
    
    def __init__(self, data_directory: str, output_file: str = "tov_mass_radius_data.csv", r_end_cm: float | None = 1.0e8):
        """
        Initialize the extractor.
        
        Args:
            data_directory: Directory containing TOV solution CSV files
            output_file: Output file for extracted M-R data
        """
        self.data_directory = data_directory
        self.output_file = output_file
        self.r_end_cm = r_end_cm
        
        # EOS file patterns (same as in stellar_plotter.py)
        self.eos_patterns = {
            'hybrid': 'tov_solution_magnetic_bps_bbp_polytrope_*.csv',
            'neutron_relativistic': 'neutron_relativistic_rhoc_*.csv',
            'electron_relativistic': 'electron_relativistic_rhoc_*.csv',
            'electron_non_relativistic': 'electron_non_relativistic_rhoc_*.csv',
            'neutron_non_relativistic': 'neutron_non_relativistic_rhoc_*.csv'
        }
        
        # Physical constants
        self.SOLAR_MASS_G = 1.989e33  # Solar mass in grams
        self.CM_TO_KM = 1e-5         # cm to km conversion
        
        # Results storage
        self.results = []
    
    def parse_cpp_encoded_density(self, density_str: str) -> float:
        """
        Parse density string using C++ encoding rules.
        
        C++ encoding converts:
        - '1.00e+16' -> '1.00pp16' (positive exponent)
        - '1.50e-03' -> '1.50p-03' (negative exponent)
        
        Args:
            density_str: Encoded density string from filename
            
        Returns:
            Central density in g/cm³
        """
        try:
            # Reverse the C++ encoding
            if 'pp' in density_str:
                # Handle positive exponent: '1.00pp16' -> '1.00e+16'
                density_str = density_str.replace('pp', 'e+')
            elif density_str.count('p') == 1:
                # Handle negative exponent: '1.50p-03' -> '1.50e-03'
                density_str = density_str.replace('p', 'e', 1)
            
            return float(density_str)
        except ValueError as e:
            logger.error(f"Failed to parse density '{density_str}': {e}")
            return None
    
    def parse_central_density_from_filename(self, filename: str, eos_type: str) -> Optional[float]:
        """
        Extract central density from filename based on EOS type.
        
        Args:
            filename: CSV filename
            eos_type: EOS type identifier
            
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
                match = re.search(rf'{eos_type}_rhoc_(.*)\.csv', base)
                if match:
                    rho_str = match.group(1)  # e.g. '5.00pp18'
                    return self.parse_cpp_encoded_density(rho_str)
                    
        except (ValueError, AttributeError) as e:
            logger.error(f"Failed to parse density from {filename}: {e}")
        
        return None
    
    def extract_mass_radius_from_csv(self, filepath: str, eos_type: str) -> Optional[Dict]:
        """
        Extract total mass and radius from a single TOV solution CSV file.
        
        This is the core extraction logic from stellar_plotter.py
        
        Args:
            filepath: Path to CSV file
            eos_type: EOS type identifier
            
        Returns:
            Dictionary with extracted data, or None if extraction fails
        """
        try:
            # Load CSV data
            df = pd.read_csv(filepath)
            
            # Check for required columns
            required_cols = ['log_m[g]', 'log_r[cm]']
            if not all(col in df.columns for col in required_cols):
                logger.error(f"Missing required columns in {filepath}")
                logger.error(f"Available columns: {list(df.columns)}")
                return None
            
            # Remove NaN values
            df = df.dropna()
            if df.empty:
                logger.error(f"No valid data in {filepath}")
                return None
            
            # Extract final (surface) values - KEY STEP!
            # The last row contains the total mass and total radius
            last_mass_log = df['log_m[g]'].iloc[-1]    # log10(total mass in grams)
            last_radius_log = df['log_r[cm]'].iloc[-1]  # log10(total radius in cm)
            
            # Convert to physical units
            mass_solar = 10**last_mass_log / self.SOLAR_MASS_G  # grams to solar masses
            radius_km = 10**last_radius_log * self.CM_TO_KM     # cm to km
            
            # Parse central density from filename
            central_density = self.parse_central_density_from_filename(filepath, eos_type)
            if central_density is None:
                logger.warning(f"Could not parse central density from {filepath}")
                central_density = np.nan

            # Skip runs that hit the outer radius instead of finding the surface
            if self.r_end_cm is not None:
                r_end_km = self.r_end_cm * self.CM_TO_KM
                # allow tiny FP wiggle; skip if we're essentially at r_end
                if np.isclose(radius_km, r_end_km, rtol=1e-3, atol=0.5):
                    logger.warning(f"{os.path.basename(filepath)} ended at r_end≈{r_end_km:.2f} km; skipping.")
                    return None

            # Verify the last row is actually the surface (bracketed crossing)
            pcol = 'log_P[dyne/cm^2]'
            if pcol in df.columns:
                last_lp = float(df[pcol].iloc[-1])
                min_lp  = float(df[pcol].min())
                has_surface = np.isclose(last_lp, min_lp, atol=5e-4)  # 5e-4 dex tolerance
                if not has_surface:
                    logger.warning(f"{os.path.basename(filepath)} did not reach surface (hit r_end); skipping.")
                    return None
            else:
                logger.warning(f"{os.path.basename(filepath)} missing '{pcol}' column; cannot verify surface reliably.")

            
            # Prepare result dictionary
            result = {
                'filename': os.path.basename(filepath),
                'eos_type': eos_type,
                'central_density_g_cm3': central_density,
                'log_central_density': math.log10(central_density) if central_density else np.nan,
                'total_mass_solar': mass_solar,
                'total_radius_km': radius_km,
                'log_mass_g': last_mass_log,
                'log_radius_cm': last_radius_log,
                'data_points': len(df)  # Number of radial points in the solution
            }
            
            logger.info(f"Extracted: {os.path.basename(filepath)} -> "
                       f"M={mass_solar:.4f} M☉, R={radius_km:.2f} km, "
                       f"ρc={central_density:.2e} g/cm³")
            
            return result
            
        except Exception as e:
            logger.error(f"Failed to process {filepath}: {e}")
            return None
    
    def discover_and_process_files(self, eos_types: Optional[List[str]] = None) -> None:
        """
        Discover and process all TOV solution files.
        
        Args:
            eos_types: List of EOS types to process. If None, processes all types.
        """
        if eos_types is None:
            eos_types = list(self.eos_patterns.keys())
        
        total_files = 0
        processed_files = 0
        
        for eos_type in eos_types:
            if eos_type not in self.eos_patterns:
                logger.warning(f"Unknown EOS type: {eos_type}")
                continue
            
            # Find files matching the pattern
            pattern = os.path.join(self.data_directory, self.eos_patterns[eos_type])
            files = glob.glob(pattern)
            
            if not files:
                logger.warning(f"No files found for {eos_type} with pattern: {pattern}")
                continue
            
            logger.info(f"Processing {len(files)} files for EOS type: {eos_type}")
            total_files += len(files)
            
            # Process each file
            for filepath in sorted(files):
                result = self.extract_mass_radius_from_csv(filepath, eos_type)
                if result:
                    self.results.append(result)
                    processed_files += 1
        
        logger.info(f"Successfully processed {processed_files}/{total_files} files")
    
    def save_results(self, sort_by: str = 'central_density_g_cm3') -> None:
        """
        Save extracted M-R data to CSV file.
        
        Args:
            sort_by: Column to sort results by
        """
        if not self.results:
            logger.error("No results to save")
            return
        
        # Convert to DataFrame
        df = pd.DataFrame(self.results)
        
        # Sort results
        if sort_by in df.columns:
            df = df.sort_values([sort_by, 'eos_type']).reset_index(drop=True)
        
        # Add some derived quantities
        df['mass_radius_ratio'] = df['total_mass_solar'] / df['total_radius_km']
        df['compactness'] = df['total_mass_solar'] * 1.477 / df['total_radius_km']  # GM/Rc²
        
        # Save to CSV
        output_path = os.path.join(os.path.dirname(self.data_directory), self.output_file)
        df.to_csv(output_path, index=False, float_format='%.6e')
        
        logger.info(f"Saved {len(df)} M-R data points to: {output_path}")
        
        # Print summary statistics
        self.print_summary_statistics(df)
    
    def print_summary_statistics(self, df: pd.DataFrame) -> None:
        """Print summary statistics for the extracted data."""
        logger.info("\n" + "="*60)
        logger.info("SUMMARY STATISTICS")
        logger.info("="*60)
        
        for eos_type in df['eos_type'].unique():
            eos_data = df[df['eos_type'] == eos_type]
            logger.info(f"\n{eos_type.replace('_', ' ').title()}:")
            logger.info(f"  Number of models: {len(eos_data)}")
            logger.info(f"  Mass range: {eos_data['total_mass_solar'].min():.4f} - "
                       f"{eos_data['total_mass_solar'].max():.4f} M☉")
            logger.info(f"  Radius range: {eos_data['total_radius_km'].min():.2f} - "
                       f"{eos_data['total_radius_km'].max():.2f} km")
            logger.info(f"  Density range: {eos_data['central_density_g_cm3'].min():.2e} - "
                       f"{eos_data['central_density_g_cm3'].max():.2e} g/cm³")
            logger.info(f"  Maximum mass: {eos_data['total_mass_solar'].max():.4f} M☉")
            logger.info(f"  Radius at max mass: "
                       f"{eos_data.loc[eos_data['total_mass_solar'].idxmax(), 'total_radius_km']:.2f} km")
        
        logger.info(f"\nOverall statistics:")
        logger.info(f"  Total models: {len(df)}")
        logger.info(f"  Mass range: {df['total_mass_solar'].min():.4f} - {df['total_mass_solar'].max():.4f} M☉")
        logger.info(f"  Radius range: {df['total_radius_km'].min():.2f} - {df['total_radius_km'].max():.2f} km")


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Extract total mass and radius from TOV solution CSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Extract M-R data from all EOS types
  python mr_extractor.py --data-dir ./data/
  
  # Extract only hybrid and neutron relativistic
  python mr_extractor.py --data-dir ./data/ --eos-types hybrid neutron_relativistic
  
  # Custom output file
  python mr_extractor.py --data-dir ./data/ --output mr_results.csv
        """
    )
    
    parser.add_argument('--data-dir', required=True, 
                       help='Directory containing TOV solution CSV files')
    parser.add_argument('--output', default='tov_mass_radius_data.csv',
                       help='Output CSV filename (default: tov_mass_radius_data.csv)')
    parser.add_argument('--eos-types', nargs='+',
                       choices=['hybrid', 'neutron_relativistic', 'electron_relativistic',
                               'electron_non_relativistic', 'neutron_non_relativistic'],
                       help='EOS types to process (default: all)')
    parser.add_argument('--sort-by', default='central_density_g_cm3',
                       choices=['central_density_g_cm3', 'total_mass_solar', 'total_radius_km'],
                       help='Column to sort results by (default: central_density_g_cm3)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--r-end-cm', type=float, default=1.0e8,
                    help='Treat solutions that end at or near this radius as non-surface (default: 1e8 cm)')
    
    args = parser.parse_args()
    
    # Configure logging level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Validate data directory
    if not os.path.exists(args.data_dir):
        logger.error(f"Data directory does not exist: {args.data_dir}")
        return
    
    # Initialize extractor
    extractor = TOVMassRadiusExtractor(args.data_dir, args.output, r_end_cm=args.r_end_cm)
    
    # Extract data
    logger.info(f"Starting M-R extraction from: {args.data_dir}")
    extractor.discover_and_process_files(args.eos_types)
    
    # Save results
    extractor.save_results(args.sort_by)
    
    logger.info("M-R extraction completed successfully!")


if __name__ == "__main__":
    main()