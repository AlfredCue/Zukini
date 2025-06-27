#!/usr/bin/env python3
"""
OptiZDRunner Module

Coordinates a multi-step optimization pipeline for NOx reduction using ZDPlasKin,
involving reaction filtering, kinetic file generation, error evaluation, and
reaction set recovery.

Author: Alfredo Cuellar Valencia
"""

import os
import shutil
import logging
import json
import subprocess
from typing import List
from pathlib import Path
from datetime import datetime

# -----------------------------
# External Module Dependencies
# -----------------------------
from src.new_kinet import run_new_kinet
from src.filter_1_threshold import run_filter_thresholds
from src.pump_input_gen import generate_pumpkin_input
from src.parse_pumpkin import run_parse_pumpkin
from src.error_calc import run_error_analysis
from src.brute_force_recovery import brute_force_recovery, extract_simulation_time
from src.sensitivity_filter import run_sensitivity_filter
from src.single_run_module import ZDPlasKinRunner


class OptiZDRunner:
    """
    Main class coordinating the ZDPlasKin optimization pipeline.
    """

    def __init__(self, config_file="global_input.txt"):
        

        # Get the project root (Zukini/)
        self.PROJECT_ROOT = Path(__file__).resolve().parents[1]
        self.HELPER_FILES = self.PROJECT_ROOT / "src/helper_files"
        self.LOG_DIR = self.PROJECT_ROOT / "logs"
        self.config_file = self.PROJECT_ROOT/"input"/config_file
        self.threshold_config_file = self.HELPER_FILES / "thresholds_config.json"
        self.active_thresholds_file = self.HELPER_FILES / "active_thresholds.json"
        self.reacciones_file = self.HELPER_FILES / "reacciones_usadas.txt"
        self.reactions_tagged_file = self.HELPER_FILES / "reactions_tagged.json"
        self.stop_reason = None
        self.recovered = False

        # Paths
        self.zdplaskin_path = self.PROJECT_ROOT / "resources/ZDPlaskin"
        self.pumpkin_path = self.PROJECT_ROOT / "resources/pumpkin-1.4"
        self.bat_files = self.PROJECT_ROOT / "src/helper_bat_scripts"

        # Parameters from config
        self.species_interest = []
        self.max_error = 0.0

        # Threshold state
        self.current_thresholds = {}
        self.thresholds_config = {}

        self.setup_logging()
        self._load_config()
        self._load_thresholds()
        self.restart_tags()

    # -----------------------------
    # Setup and Config Loading
    # -----------------------------

    def setup_logging(self):
        """
        Configure logging to both file and stdout.
        """
        log_file = self.LOG_DIR / f"{Path(__file__).stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler(log_file, mode='w',encoding='utf-8'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def _load_config(self):
        """
        Load simulation parameters from the config file.
        Expected keys:
            - interest: NO, NO2,...
            - max_acum_erro: <float>
        """
        try:
            with open(self.config_file, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            raise RuntimeError(f"Failed to read config file '{self.config_file}': {e}")

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("interest:"):
                self.species_interest = [s.strip() for s in line.split(":", 1)[1].split(",")]
            elif line.startswith("max_acum_erro:"):
                self.max_error = float(line.split(":", 1)[1].strip())

        if not self.species_interest or self.max_error == 0.0:
            raise ValueError("Missing or invalid parameters in configuration file.")

    def _load_thresholds(self):
        """
        Load initial and step values for each threshold group.
        """
        try:
            with open(self.threshold_config_file, 'r') as f:
                self.thresholds_config = json.load(f)
            self.current_thresholds = {
                group: data["initial"] for group, data in self.thresholds_config.items()
            }
        except Exception as e:
            raise RuntimeError(f"Failed to load thresholds: {e}")

    def restart_tags(self):
        """
        Reset all reaction tags to False.
        """
        try:
            with open(self.reactions_tagged_file, 'r') as f:
                data = json.load(f)
            for key in data:
                data[key] = False
            with open(self.reactions_tagged_file, 'w') as f:
                json.dump(data, f, indent=4)
        except Exception as e:
            self.logger.error("Failed to reset reaction tags: %s", e)

    # -----------------------------
    # Support Methods
    # -----------------------------

    def run_generation(self, species_list):
        """
        Expand species tree using Pumpkin and parsing modules.
        """
        QT1_FILES = self.PROJECT_ROOT / "Data/QtFiles/Qt_1"
        new_species = set()
        for species in species_list:
            generate_pumpkin_input(species,QT1_FILES / "qt_species_list.txt",self.pumpkin_path / "src/Examples/Sim")
            self.run_script([os.path.join(self.pumpkin_path, "run_pumpkin.bat")], cwd=self.pumpkin_path)
            run_parse_pumpkin(QT1_FILES / "qt_reactions_list.txt", self.pumpkin_path / "pumpkin_output.txt",self.HELPER_FILES / "reactions_tagged.json",self.HELPER_FILES / "species_to_expand.json")

            EXPAND = self.HELPER_FILES / "species_to_expand.json"

            if os.path.exists(EXPAND):
                with open(EXPAND, "r") as f:
                    generated = json.load(f)
                    new_species.update(generated)

            try:
                os.remove(self.pumpkin_path / "pumpkin_output.txt")
            except FileNotFoundError:
                pass  

        return new_species

    def update_qtfiles(self, src_root: Path, pumpkin_root: Path, logger: logging.Logger):
        """
        Copy QtFiles from src_root to the Pumpkin simulation folder.
        """
        dest = pumpkin_root / "src" / "Examples" / "Sim"
        logger.info("[QtFiles] Refreshing Qt files %s", dest)
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src_root, dest, dirs_exist_ok=True)
        logger.info("[QtFiles] Qt files updated successfully.")

    def check_errors(self, path: str, threshold: float) -> bool:
        """
        Check if peak error exceeds threshold in a given log file.
        """
        try:
            with open(path, 'r') as f:
                for line in f:
                    if "peak:" in line:
                        peak = float(line.split("peak:")[1].split("%")[0].strip())
                        if peak > threshold:
                            self.logger.info(f"Exceeded error threshold: {peak}% > {threshold}%")
                            return True
        except Exception as e:
            raise RuntimeError(f"Failed to parse peak errors: {e}")
        return False

    def count_reactions(self, path: str) -> int:
        """
        Count number of active reactions in a reaction file.
        """
        try:
            with open(path, 'r') as f:
                return sum(1 for line in f if line.strip() and not line.startswith("#"))
        except Exception:
            return -1

    def run_script(self, command: List[str], cwd: str = None):
        """
        Run an external command as subprocess.
        """
        try:
            self.logger.info(f"Running script: {' '.join(command)}")
            subprocess.run(command, cwd=cwd, shell=False, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Script failed with error: {e}")
            raise
        except FileNotFoundError:
            self.logger.error(f"Script not found: {' '.join(command)}")
            raise

    # -----------------------------
    # Main Optimization Loop
    # -----------------------------

    #This is the main method that operates the optimization process

    def optimize(self, compiler="gfortran"):
        """
        Main optimization loop:
        - Expands species tree.
        - Applies filtering thresholds.
        - Runs ZDPlasKin simulations.
        - Validates results and performs recovery if needed.
        - Iteratively increases thresholds until stopping criteria.
        """

        #BEGIN OF PHASE 1 (Pick up)
        current_species = set(self.species_interest)
        processed_species = set()
        generations = 5

        #Base run for reference
        runner = ZDPlasKinRunner(compiler=compiler)
        runner.run()

        old_time = extract_simulation_time(self.PROJECT_ROOT / "resources/ZDPlaskin/simulation_output.txt")
        recovered_output_file = self.HELPER_FILES / "recovered_reactions.json"
        
        QT1_FILES = self.PROJECT_ROOT / "Data/QtFiles/Qt_1"
        QT2_FILES = self.PROJECT_ROOT / "Data/QtFiles/Qt_2"

        self.update_qtfiles(
            src_root=QT1_FILES,
            pumpkin_root= self.pumpkin_path, 
            logger=self.logger
        )

        #Cycle for Pumpkin chemical reaction network generation
        for gen in range(generations):
            self.logger.info(f"Generation {gen} - processing {len(current_species)} species")
            current_species -= processed_species
            if not current_species:
                self.logger.info("No further species to expand tree")
                break
            new_species = self.run_generation(current_species)
            processed_species.update(current_species)
            current_species = new_species

        #Sensitivity filter
        run_sensitivity_filter(self.species_interest)

        #END OF PHASE 1

        previous_count = -1
        counter = 1
        stop = False
        self.recovered = False
        i = 0

        with open(recovered_output_file, 'w') as f:
            json.dump({}, f)

        # Qt_2 -> Qt_1 baseline
        for file_name in os.listdir(QT1_FILES):
            shutil.copy2(
                os.path.join(QT1_FILES, file_name),
                os.path.join(QT2_FILES, file_name)
            )
        self.logger.info("Initial Qt files copied from Qt_2 to Qt_1 for baseline reference.")

        #While loop for Elimination process (PHASE 2)

        while not stop:
            try:
                with open(self.active_thresholds_file, 'w') as f:
                    json.dump(self.current_thresholds, f, indent=2)

                # Backup Qt_2
                backup_dir = self.PROJECT_ROOT / "Data/temp/qt_backup"
                backup_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["qt_rates.txt", "qt_species_list.txt", "qt_reactions_list.txt"]:
                    shutil.copy2(QT2_FILES / fname, backup_dir / fname)

                filtered, all_protected = run_filter_thresholds(self.active_thresholds_file)
                run_new_kinet()
                self.run_script([os.path.join(self.bat_files, "move_kinet.bat")])

                count = self.count_reactions(self.reacciones_file)
                self.logger.info(f"Reactions in current iteration: {count}")

                if filtered:
                    runner.run()
                    run_error_analysis(self.species_interest)
                    error_exceeded = self.check_errors(self.HELPER_FILES / "errores_acumulados.txt", self.max_error)

                    self.recovered = brute_force_recovery(
                        reference_time_sec=old_time,
                        reduced_log= self.zdplaskin_path / "simulation_output.txt",
                        filtered_reactions_file= self.HELPER_FILES / "reacciones_filtradas.txt",
                        reactions_tagged_file=self.reactions_tagged_file,
                        recovered_output_file=recovered_output_file,
                        force=error_exceeded
                    )
                    old_time = extract_simulation_time(self.zdplaskin_path / "simulation_output.txt")

                if self.recovered:
                    self.logger.info("Brute-force recovery succeeded. Restoring Qt_2 and re-filtering...")
                    for fname in ["qt_rates.txt", "qt_species_list.txt", "qt_reactions_list.txt"]:
                        shutil.copy2(backup_dir / fname, QT2_FILES / fname)
                else:
                    self.logger.info("No recovery needed or recovery failed.")
                    if all_protected:
                        stop = True
                        self.logger.info("No more reactions can be eliminated, stopping optimization....")
                    else:
                        shutil.copy(self.zdplaskin_path / "kinet.inp", self.PROJECT_ROOT / "output/kinet-final.inp")

                # Increase thresholds for next iteration
                for group in self.current_thresholds:
                    step = self.thresholds_config[group]["step"]
                    mult = 1 if count != previous_count else 100
                    self.current_thresholds[group] += step * mult * counter

                previous_count = count
                counter = counter ** 1.1 + 1
                i += 1
                self.logger.info("Threshold increased, iteration number: %d", i)

            except Exception as e:
                self.logger.critical(f"Unrecoverable error: {e}")
                self.stop_reason = "exception"
                break

        self.logger.info("Optimization loop ended.")
