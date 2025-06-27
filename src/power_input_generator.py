#!/usr/bin/env python3
"""
Power Input Generator

This module interfaces with MATLAB to generate power input files (`data_in.dat` and `run_info.txt`)
for a specified duty cycle. It writes output into a configurable directory (`power_info/` by default).
"""

import subprocess
from pathlib import Path
import logging

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)

# ----------------------------------------------------------------------
# Power Input Generator Class
# ----------------------------------------------------------------------
class PowerInputGenerator:
    """
    Manages the generation of power input files using an external MATLAB script.

    The script must accept two arguments:
    1. duty_cycle (float)
    2. output directory (string path)

    Example MATLAB call constructed:
        main(0.3, 'power_info');exit;
    """

    def __init__(
        self,
        matlab_main_fn: str = "main",
        output_dir: Path = Path("power_info"),
        matlab_path: str = "matlab"
    ):
        """
        Args:
            matlab_main_fn (str): Name of the MATLAB function to call.
            output_dir (Path): Directory to store output files.
            matlab_path (str): Path to MATLAB executable (e.g., with batch mode support).
        """
        PROJECT_ROOT = Path(__file__).resolve().parents[1]
        self.matlab_main_fn = matlab_main_fn
        self.output_dir = PROJECT_ROOT / "temp/power_info"
        self.matlab_path = matlab_path
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, duty_cycle: float) -> Path:
        """
        Executes the MATLAB script to generate input for the given duty cycle.

        Args:
            duty_cycle (float): Value between 0 and 1.

        Returns:
            Path: Folder containing data_in.dat and run_info.txt.

        Raises:
            ValueError: If duty cycle is outside (0, 1].
            RuntimeError: If MATLAB returns an error.
            FileNotFoundError: If expected output files are missing.
        """
        if not (0 < duty_cycle <= 1):
            raise ValueError("Duty cycle must be between 0 and 1.")

        # Clean output folder if needed
        if self.output_dir.exists():
            for file in self.output_dir.glob("*"):
                file.unlink()
        else:
            self.output_dir.mkdir(parents=True)

        matlab_cmd = f"{self.matlab_main_fn}({duty_cycle:.12g}, '{self.output_dir.as_posix()}');exit;"
        LOGGER.info("Running MATLAB for D = %.4f", duty_cycle)

        result = subprocess.run(
            [self.matlab_path, "-batch", matlab_cmd],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            LOGGER.error("MATLAB stdout: %s", result.stdout)
            LOGGER.error("MATLAB stderr: %s", result.stderr)
            raise RuntimeError(f"MATLAB run failed for D = {duty_cycle:.4f}")

        data_file = self.output_dir / "data_in.dat"
        info_file = self.output_dir / "run_info.txt"

        if not data_file.exists() or not info_file.exists():
            raise FileNotFoundError("Expected output files not found in power_info/")

        LOGGER.info("MATLAB output ready in %s", self.output_dir)
        return self.output_dir

# ----------------------------------------------------------------------
# CLI Testing (optional)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    gen = PowerInputGenerator()
    try:
        gen.run(0.3)  # Example test duty cycle
    except Exception as e:
        LOGGER.exception("Failed to generate power input: %s", e)
