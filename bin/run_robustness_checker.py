#!/usr/bin/env python3
"""
Robustness Analysis Orchestrator

Performs a full factorial study across:
- Temperatures (T)
- Gas mixtures (NO, O2 ppm)
- Duty cycles (D)

For each combination, runs ZDPlasKin with:
- Original kinetic file
- Filtered set 1
- Filtered set 2

Results are stored under `robustness_analysis_results/` and summarized into a final CSV.
"""

import itertools
import json
import logging
import shutil
from pathlib import Path
from datetime import datetime

from src.gas_mix_builder import GasMixBuilder
from src.const_temp_builder import TemperatureBuilder
from src.single_run_module import ZDPlasKinRunner
from src.plot_densities_robustness import ZDPlasKinAnalyzer
from src.power_input_generator import PowerInputGenerator
from src.brute_force_recovery import extract_simulation_time
from src.robustness_analysis_module import RobustnessResultsAggregator

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = PROJECT_ROOT / "logs"
log_file = LOG_DIR / f"{Path(__file__).stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode='w',encoding='utf-8')]
)
LOGGER = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# Parameters and Paths
# ----------------------------------------------------------------------

#This values can be changed in case certain conditions want to be checked
TEMP_VALUES = [300, 400, 500]
DUTY_VALUES = [0.1, 0.5, 0.9]
GAS_VALUES = [
    {"NO": 800, "O2": 25000},
    {"NO": 900, "O2": 30000},
    {"NO": 1000, "O2": 35000},
]
KINETIC_SETS = [
    (PROJECT_ROOT / "data/Reduced_Sets/kinet_original.inp", "original"),
    (PROJECT_ROOT / "data/Reduced_Sets/kinet_filtered1.inp", "filtered1"),
    (PROJECT_ROOT / "data/Reduced_Sets/kinet_filtered2.inp", "filtered2"),
]

FORTRAN_INPUT = PROJECT_ROOT / "resources/ZDPlaskin/dbd0D.F90"
TEMP_JSON = PROJECT_ROOT / "src/helper_files/const_temp_param.json"
GAS_JSON = PROJECT_ROOT / "src/helper_files/gas_mixture.json"

ZDPLASKIN_DIR = PROJECT_ROOT / "resources/ZDPlaskin"
QT1_DIR = PROJECT_ROOT / "data/QtFiles/Qt_1"
QT2_DIR = PROJECT_ROOT / "data/QtFiles/Qt_2"
RESULTS_DIR = PROJECT_ROOT / "output/robustness_analysis_results"

# ----------------------------------------------------------------------
# Utility: Safe move (with parent creation)
# ----------------------------------------------------------------------
def safe_move(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.move(str(src), str(dst))
    else:
        LOGGER.warning(f"File or directory not found for move: {src}")

# ----------------------------------------------------------------------
# Main Orchestration Function
# ----------------------------------------------------------------------
def run_parametric_study():
    """
    Performs the full robustness check pipeline.
    """
    runner = ZDPlasKinRunner(compiler="gfortran")
    analyzer = ZDPlasKinAnalyzer()
    power_generator = PowerInputGenerator(matlab_path="C:/Program Files/MATLAB/R2024a/bin/matlab")

    for temp, gas, duty in itertools.product(TEMP_VALUES, GAS_VALUES, DUTY_VALUES):
        try:
            # Generate folder name
            gas_str = f"G_NO_{gas['NO']}_O2_{gas['O2']}"
            duty_str = f"Dc_{str(duty).replace('.', '_')}"
            config_name = f"T{temp}_{gas_str}_{duty_str}"
            output_base = RESULTS_DIR / config_name

            LOGGER.info("Starting configuration: %s", config_name)

            # Update temperature (via JSON + builder)
            TEMP_JSON.write_text(json.dumps({"gas_temperature": temp}))
            TemperatureBuilder(TEMP_JSON, FORTRAN_INPUT, FORTRAN_INPUT).update()

            # Update gas mixture
            GAS_JSON.write_text(json.dumps(gas))
            GasMixBuilder(GAS_JSON, FORTRAN_INPUT, FORTRAN_INPUT).update()

            # Generate and copy power input
            generated_path = power_generator.run(duty_cycle=duty)
            shutil.move(str(generated_path), ZDPLASKIN_DIR / "data_in.dat")

            # Loop over kinetic sets
            for kinetic_file, label in KINETIC_SETS:
                shutil.copyfile(kinetic_file, ZDPLASKIN_DIR / "kinet.inp")
                runner.run()

                if label == "original":
                    # Move Qt_2 to Qt_1 as reference
                    if QT1_DIR.exists():
                        shutil.rmtree(QT1_DIR)
                    shutil.move(str(QT2_DIR), str(QT1_DIR))
                    sim_time_full = extract_simulation_time(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                else:
                    sim_time_reduced = extract_simulation_time(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                    result_path = output_base / label

                    score = analyzer.analyze(
                        qt1_dir= QT1_DIR,
                        qt2_dir= QT2_DIR,
                        priority_species=["NO", "NO2"],
                        sim_time_full=sim_time_full,
                        sim_time_reduced=sim_time_reduced
                    )
                    

                    # Save results
                    if (img_dir := PROJECT_ROOT / "images").exists():
                        shutil.copytree(img_dir, result_path, dirs_exist_ok=True)
                        shutil.move(RESULTS_DIR / "metrics_summary.txt", RESULTS_DIR / config_name / label)
                    LOGGER.info("Trust score for %s: %.2f", label, score)

        except Exception as e:
            LOGGER.error("Failed configuration %s due to: %s", config_name, str(e), exc_info=True)

    # Compile all results into CSV
    agg = RobustnessResultsAggregator()
    df = agg.compile(
        root=RESULTS_DIR,
        out_csv= RESULTS_DIR / "compiled_results.csv"
    )

# ----------------------------------------------------------------------
# Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    run_parametric_study()
