#!/usr/bin/env python3
"""
ZDPlasKin Runner Module

Provides a unified interface for launching ZDPlasKin simulations using different compilers.
After running, it also triggers Qt output file relocation via move_qtfile.bat.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Literal

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)

# ----------------------------------------------------------------------
# ZDPlasKin Execution Wrapper
# ----------------------------------------------------------------------
class ZDPlasKinRunner:
    """
    Main class to run a ZDPlasKin simulation using the appropriate batch script,
    depending on the selected compiler ("gfortran" or "ifort").
    """

    def __init__(self, compiler: Literal["gfortran", "ifort"]):
        """
        Args:
            compiler (str): Compiler used for build ("gfortran" or "ifort").
        """
        self.compiler = compiler
        self.PROJECT_ROOT = Path(__file__).resolve().parents[1]
        self.zdplaskin_path = self.PROJECT_ROOT / "resources/ZDPlaskin"

    def _run_script(self, command: List[str], cwd: str = None):
        """
        Run an external batch script as a subprocess.

        Args:
            command (List[str]): List with command path and arguments.
            cwd (str): Working directory (optional).

        Raises:
            CalledProcessError or FileNotFoundError if script fails.
        """
        try:
            LOGGER.info(f"Running script: {' '.join(command)}")
            subprocess.run(command, cwd=cwd, shell=False, check=True)
        except subprocess.CalledProcessError as e:
            LOGGER.error(f"Script failed with error: {e}")
            raise
        except FileNotFoundError:
            LOGGER.error(f"Script not found: {' '.join(command)}")
            raise

    def run(self):
        """
        Run ZDPlasKin simulation and move Qt files afterward.
        Chooses batch script based on selected compiler.
        """
        if self.compiler == "ifort":
            script = os.path.join(self.zdplaskin_path, "run_zdplaskin.bat")
        elif self.compiler == "gfortran":
            script = os.path.join(self.zdplaskin_path, "run_zdplaskin_gfort.bat")
        else:
            raise ValueError(f"Unsupported compiler: {self.compiler}")

        self._run_script([script], cwd=self.zdplaskin_path)
        self._run_script([os.path.join(self.PROJECT_ROOT, "src/helper_bat_scripts/move_qtfile.bat")])
        LOGGER.info("Simulation run completed successfully.")

# ----------------------------------------------------------------------
# CLI Usage Example (Debug/Test)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    runner = ZDPlasKinRunner(compiler="gfortran")
    runner.run()
