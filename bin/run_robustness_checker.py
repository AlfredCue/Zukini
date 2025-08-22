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
import time
from pathlib import Path
from datetime import datetime

from src.gas_mix_builder import GasMixBuilder
from src.const_temp_builder import TemperatureBuilder
from src.single_run_module import ZDPlasKinRunner
from src.plot_densities_robustness import ZDPlasKinAnalyzer
from src.power_input_generator import PowerInputGenerator
from src.brute_force_recovery import extract_simulation_time
from src.brute_force_recovery import extract_selectivity
from src.robustness_analysis_module import RobustnessResultsAggregator
from bin.run_optizd import clear_qt_directories
from src.plot_densities_robustness import extract_last_value, load_species_names, load_densities


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
TEMP_VALUES = [300,350,400,450,500,550,600]
DUTY_VALUES = [0.5]
GAS_VALUES = [
    {"NO": 800, "O2": 0},
    {"NO": 800, "O2": 25000}
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
    power_generator = PowerInputGenerator(matlab_path="C:/Program Files/MATLAB/R2025a/bin/matlab")

    qt_dirs = [
        PROJECT_ROOT / "Data/QtFiles/Qt_1",
        PROJECT_ROOT / "Data/QtFiles/Qt_2"
    ]
    clear_qt_directories(qt_dirs)

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
            data_in_path = ZDPLASKIN_DIR / "data_in.dat"
            if data_in_path.exists():
                data_in_path.unlink()
            shutil.move(generated_path / "data_in.dat", ZDPLASKIN_DIR)

            if QT1_DIR.exists():
                shutil.rmtree(QT1_DIR)
            # Loop over kinetic sets
            for kinetic_file, label in KINETIC_SETS:
                shutil.copyfile(kinetic_file, ZDPLASKIN_DIR / "kinet.inp")
                runner.run()
                if label == "original":
                    t_full, dens_full = load_densities(QT1_DIR / "qt_densities.txt")
                    species_names = load_species_names(QT1_DIR / "qt_species_list.txt")
                    sim_time_full = extract_simulation_time(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                    selectivity_full = extract_selectivity(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                    o3_full = extract_last_value("O3",species_names=species_names,density_matrix=dens_full)
                else:
                    t_red, dens_red = load_densities(QT2_DIR / "qt_densities.txt")
                    sim_time_reduced = extract_simulation_time(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                    selectivity_reduced = extract_selectivity(PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
                    o3_reduced = extract_last_value("O3",species_names=species_names,density_matrix=dens_red)
                    result_path = output_base / label

                    score = analyzer.analyze(
                        qt_full= QT1_DIR,
                        qt_reduced= QT2_DIR,
                        priority_species=["NO", "NO2"],
                        sim_time_full=sim_time_full,
                        sim_time_reduced=sim_time_reduced,
                        selectivity_full=selectivity_full,
                        selectivity_reduced=selectivity_reduced,
                        o3_full=o3_full,
                        o3_reduced=o3_reduced
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
