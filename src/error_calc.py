#!/usr/bin/env python3
"""
Error Analysis Module

Compares density data for a given set of species between two QtPlasKin outputs
(Qt_1 and Qt_2), and computes:
- Mean relative error
- Relative RMSE (normalized by max absolute value of original)

Results are saved to errores_acumulados.txt and printed to stdout.
"""

import numpy as np
import sys
import os
import logging
from typing import List, Tuple
import argparse
from pathlib import Path
from datetime import datetime

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = PROJECT_ROOT / "logs"
log_file = LOG_DIR / f"{Path(__file__).stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_file, mode='w',encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# Data Loading Functions
# ----------------------------------------------------------------------
def load_densities(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load time and density values from qt_densities.txt.

    Returns:
        times (np.ndarray), densities (np.ndarray [n_timesteps, n_species])
    """
    try:
        data = np.loadtxt(file_path, skiprows=1)
        return data[:, 0], data[:, 1:]
    except Exception as e:
        raise RuntimeError(f"Failed to load density file '{file_path}': {e}")

def load_species_names(file_path: str) -> List[str]:
    """
    Load species names from qt_species_list.txt (ignores indices).

    Returns:
        List of species names in order.
    """
    try:
        with open(file_path, 'r') as f:
            return [line.strip().split(maxsplit=1)[1] for line in f if line.strip()]
    except Exception as e:
        raise RuntimeError(f"Failed to load species list '{file_path}': {e}")

# ----------------------------------------------------------------------
# Extraction and Metrics
# ----------------------------------------------------------------------
def extract_density_by_name(name: str, all_names: List[str], all_densities: np.ndarray) -> np.ndarray:
    """
    Extract density time series for a given species name.

    Args:
        name: species to look for.
        all_names: list of species.
        all_densities: density matrix (time x species)

    Returns:
        density time series (np.ndarray)
    """
    if name not in all_names:
        raise ValueError(f"Species '{name}' not found in species list.")
    index = all_names.index(name)
    return all_densities[:, index]

def compute_relative_errors(dens1: np.ndarray, dens2: np.ndarray) -> Tuple[float, float]:
    """
    Compute average relative error and relative RMSE (normalized by max of dens1).

    Args:
        dens1: reference density
        dens2: reduced set density

    Returns:
        avg_error: mean relative error
        relative_rmse: RMSE / max(abs(dens1))
    """
    epsilon = 1e-20 * np.max(np.abs(dens1))  # avoid divide-by-zero
    relative_errors = np.abs(dens1 - dens2) / (np.abs(dens1) + epsilon)

    avg_error = np.mean(relative_errors)

    rmse = np.sqrt(np.mean((dens1 - dens2) ** 2))
    norm = np.max(np.abs(dens1)) + epsilon
    relative_rmse = rmse / norm

    return avg_error, relative_rmse

# ----------------------------------------------------------------------
# Main Entry Function
# ----------------------------------------------------------------------
def run_error_analysis(target_species: List[str]) -> None:
    """
    Main routine for comparing densities between Qt_1 and Qt_2 for a list of species.
    Outputs results to errores_acumulados.txt.
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    HELPER_FILES = PROJECT_ROOT / "src/helper_files"
    output_file = HELPER_FILES / "errores_acumulados.txt"

    # File paths
    QT1_FILES = PROJECT_ROOT / "Data/QtFiles/Qt_1"
    QT2_FILES = PROJECT_ROOT / "Data/QtFiles/Qt_2"
    path_dens1 = QT1_FILES / "qt_densities.txt"
    path_names1 = QT1_FILES / "qt_species_list.txt"
    path_dens2 = QT2_FILES / "qt_densities.txt"
    path_names2 = QT2_FILES / "qt_species_list.txt"

    try:
        logger.info("Loading density data...")
        times_1, densities_1 = load_densities(path_dens1)
        names_1 = load_species_names(path_names1)

        times_2, densities_2 = load_densities(path_dens2)
        names_2 = load_species_names(path_names2)

        if not np.allclose(times_1, times_2, rtol=1e-1, atol=1e-1):
            raise ValueError("Time vectors from the two density files do not match.")

    except Exception as e:
        logger.critical(f"Fatal error during input loading: {e}")
        raise

    # Compare and output errors
    with open(output_file, "w") as f_out:
        for species in target_species:
            try:
                dens1 = extract_density_by_name(species, names_1, densities_1)
                dens2 = extract_density_by_name(species, names_2, densities_2)

                avg_err, rel_rmse = compute_relative_errors(dens1, dens2)
                avg_pct = avg_err * 100
                rmse_pct = rel_rmse * 100

                logger.info(f"{species}: Avg Error = {avg_pct:.2f}%, Relative RMSE = {rmse_pct:.2f}%")
                f_out.write(f"{species} = {avg_pct:.2f}% (Relative RMSE: {rmse_pct:.2f}%)\n")

            except ValueError as ve:
                logger.warning(ve)
                f_out.write(f"{species} = Species not found.\n")
            except Exception as e:
                logger.error(f"Error processing {species}: {e}")
                f_out.write(f"{species} = Error during computation.\n")

    logger.info(f"Error computation completed. Output written to {output_file}")

# ----------------------------------------------------------------------
# CLI Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute relative errors for selected species from Qt_1 vs Qt_2 density data."
    )
    parser.add_argument(
        "-s", "--species", nargs="+", required=True,
        help="List of species to analyze (e.g., -s NO NO2 O N)"
    )
    args = parser.parse_args()
    run_error_analysis(args.species)
