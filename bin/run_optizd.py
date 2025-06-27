#!/usr/bin/env python3
"""
Run Optimization with OptiZD

This script:
1. Cleans previous Qt simulation output directories.
2. Copies the original kinetic file into the working directory.
3. Launches the optimization pipeline using OptiZDRunner.

Author: Alfredo Cuellar Valencia
"""

import shutil
import os
import logging
from src.opti_orchestrator import OptiZDRunner
from pathlib import Path
from datetime import datetime


# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = PROJECT_ROOT / "logs"
log_file = LOG_DIR / f"{Path(__file__).stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode='w',encoding='utf-8')]
)

# ----------------------------------------------------------------------
# Utility Function to Clean Qt Directories
# ----------------------------------------------------------------------
def clear_qt_directories(directories: list[str]) -> None:
    """
    Delete all contents inside the specified Qt output directories.

    Args:
        directories (list[str]): List of directory paths to clear.
    """
    for dir_path in directories:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            file_path = os.path.join(dir_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                LOGGER.info(f"Deleted: {file_path}")
            except Exception as e:
                LOGGER.warning(f"Failed to delete {file_path}: {e}")

# ----------------------------------------------------------------------
# Main Pipeline Function
# ----------------------------------------------------------------------
def run_optizd_pipeline(
    source_path: str = "Data/Reduced_sets/kinet_original.inp",
    dest_path: str = "ZDPlaskin/kinet.inp",
    config_file: str = "global_input.txt",
    compiler: str = "gfortran"
) -> None:
    """
    Run the full OptiZD optimization pipeline.

    Args:
        source_path (str): Path to the original kinetic file.
        dest_path (str): Destination path where the kinetic file is copied.
        config_file (str): Path to the configuration input file.
        compiler (str): Compiler to use during optimization ("gfortran", etc.)
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    # Step 1: Clean previous Qt outputs
    qt_dirs = [
        PROJECT_ROOT / "Data/QtFiles/Qt_1",
        PROJECT_ROOT / "Data/QtFiles/Qt_2"
    ]
    LOGGER.info("Cleaning Qt directories...")
    clear_qt_directories(qt_dirs)

    # Step 2: Copy and rename kinetic input file
    if not os.path.exists(source_path):
        raise FileNotFoundError(f"Source file not found: {source_path}")

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    shutil.copyfile(source_path, dest_path)
    LOGGER.info("Copied %s to %s", source_path, dest_path)

    # Step 3: Run optimization using OptiZDRunner
    LOGGER.info("--------INITIATING OPTIMIZATION ALGORITHM--------")
    runner = OptiZDRunner(config_file=config_file)
    runner.optimize(compiler=compiler)
    LOGGER.info("Optimization complete using %s compiler.", compiler)
    LOGGER.info("--------OPTIMIZATION ENDED SUCCESFULLY--------")

# ----------------------------------------------------------------------
# Script Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    
    run_optizd_pipeline(source_path= PROJECT_ROOT / "Data/Reduced_sets/kinet_original.inp",
                        dest_path= PROJECT_ROOT / "resources/ZDPlaskin/kinet.inp",
                        config_file= PROJECT_ROOT / "input/global_input.txt",
                        compiler= "gfortran") #For intel is "ifort"
