#!/usr/bin/env python3
"""
Sweep Orchestrator

Coordinates a full pipeline for each simulation configuration:
1. Power input generation via MATLAB
2. ZDPlasKin simulation
3. Density analysis
4. Results archival into structured experiment folders
"""

from pathlib import Path
import shutil
import logging
import re
from datetime import datetime

from src.power_input_generator import PowerInputGenerator
from src.density_analysis_module import DensityAnalyzer
from src.single_run_module import ZDPlasKinRunner
from src.brute_force_recovery import extract_simulation_time

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
# Constants
# ----------------------------------------------------------------------
EXPERIMENT_ROOT = PROJECT_ROOT / "output/experiments"
QT_SOURCE_DIR = PROJECT_ROOT / "data/QtFiles/Qt_2"

# ----------------------------------------------------------------------
# Utility: Parse run_info.txt for metadata
# ----------------------------------------------------------------------
def parse_run_info(info_path: Path) -> dict:
    """
    Extract simulation metadata from run_info.txt.

    Returns:
        dict with keys: D (duty), I (current), F (Hz), P (power)
    """
    with open(info_path, 'r') as f:
        lines = f.readlines()

    data = {}
    for line in lines:
        if "Duty cycle" in line:
            data['D'] = float(re.search(r"[\d.]+", line).group())
        elif "Current" in line:
            data['I'] = float(re.search(r"[\d.]+", line).group())
        elif "Frequency" in line:
            data['F'] = int(re.search(r"\d+", line).group())
        elif "Power" in line:
            data['P'] = float(re.search(r"[\d.]+", line).group())

    return data

# ----------------------------------------------------------------------
# Utility: Create structured output folder
# ----------------------------------------------------------------------
def build_experiment_folder(base_dir: Path, info: dict) -> Path:
    """
    Build a unique folder name based on simulation parameters.

    Example:
        Duty__0_1__F_80000HZ__I_2A__P_20W
    """
    folder_name = (
        f"Duty__{str(info['D']).replace('.', '_')}"
        f"__F_{info['F']}HZ"
        f"__I_{int(round(info['I']))}A"
        f"__P_{int(round(info['P']))}W"
    )
    full_path = base_dir / folder_name
    full_path.mkdir(parents=True, exist_ok=True)
    return full_path

# ----------------------------------------------------------------------
# Main Pipeline
# ----------------------------------------------------------------------
def run_simulation_pipeline(duty_cycle: float, ts: float, qt_output_dir: Path, experiment_id: int) -> float:
    """
    Full run pipeline for a given duty cycle and time step.

    Args:
        duty_cycle (float): Duty cycle (e.g. 0.35).
        ts (float): Switching period (e.g. 0.0000125).
        qt_output_dir (Path): Directory containing output QtFiles.
        experiment_id (int): Experiment identifier (used in logs only).

    Returns:
        float: Simulation efficiency (from density analysis).
    """
    LOGGER.info("Starting simulation pipeline for duty cycle %.2f%%", duty_cycle * 100)

    #If wanted the user can skeep the step one and offer a folder with a series of input "data_in" files, only 
    #requirement is to move the data_in file to the working directory of ZDPlaskin the rest of the analyzis can be handled

    # Step 1: Generate power input using MATLAB
    # Depending on the computer the direction of the Matlab would have to be mofified
    # For the script to run the matlab requires the conection modules to PLECS and Simulink
    #----------------------------------------------------------------------------------------------
    #This is the location of the matlab executable, it has to be ajusted according to the users computer, it doesn't fetch for it 
    #automatically
    generator = PowerInputGenerator(matlab_path="C:/Program Files/MATLAB/R2025a/bin/matlab.exe")
    #----------------------------------------------------------------------------------------------
    power_dir = generator.run(duty_cycle)
    run_info_path = power_dir / "run_info.txt"
    data_in_path = power_dir / "data_in.dat"

    # Step 2: Parse simulation metadata and create experiment folder
    info = parse_run_info(run_info_path)
    output_dir = build_experiment_folder(EXPERIMENT_ROOT, info)

    # Step 3: Run ZDPlasKin simulation
    shutil.copy(data_in_path, "ZDPlaskin")  # move data_in.dat into working directory
    runner = ZDPlasKinRunner(compiler="gfortran")
    runner.run()

    # Step 4: Post-process with density analysis
    sim_time = extract_simulation_time("ZDPlaskin/simulation_output.txt")
    analyzer = DensityAnalyzer(qt_dir=qt_output_dir, ts=ts, output_dir=output_dir, sim_time=sim_time)
    efficiency = analyzer.run()

    # Step 5: Save all relevant output to the experiment folder
    shutil.copy2(run_info_path, output_dir / "run_info.txt")
    shutil.copy2(data_in_path, output_dir / "data_in.dat")

    for qtfile in qt_output_dir.glob("qt_*.txt"):
        shutil.copy2(qtfile, output_dir / qtfile.name)

    LOGGER.info("Results stored in: %s", output_dir)
    LOGGER.info("Simulation efficiency: %.6f", efficiency)

    return efficiency

# ----------------------------------------------------------------------
# CLI Entry Point for Testing
# ----------------------------------------------------------------------
if __name__ == "__main__":
    duty_cycle_fraction = 0.35
    switching_period = 5e-6
    qt_dir_path = Path("Data/QtFiles/Qt_2")

    run_simulation_pipeline(
        duty_cycle=duty_cycle_fraction,
        ts=switching_period,
        qt_output_dir=qt_dir_path,
        experiment_id=1
    )
