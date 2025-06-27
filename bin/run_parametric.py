#!/usr/bin/env python3
"""
Batch Parametric Runner

Sweeps over a range of duty cycles and runs full ZDPlasKin simulations
using `run_simulation_pipeline()` from `sweep_orchestrator`.

Designed for systematic, high-resolution parametric studies.
"""

import numpy as np
from pathlib import Path
import logging
from src.sweep_orchestrator import run_simulation_pipeline
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
# Main Execution Loop
# ----------------------------------------------------------------------
def main():
    """
    Main function to loop through a series of duty cycle values,
    invoking the pipeline for each value.
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    #This fixed switch period has to be manually modified according to the frecuency, and its asummed the 
    #power input is constant, meaning no burst mode
    ts = 0.0000125  # Fixed switching period (80 kHz) VALUE CAN CHANGE
    qt_output_dir = PROJECT_ROOT / "data/QtFiles/Qt_2" # Target directory for Qt output files
    duty_values = np.arange(0.1, 0.9, 0.2)     # Duty cycles: (initial, final, step)

    for idx, dc in enumerate(duty_values, start=1):
        LOGGER.info("Running experiment %d for duty cycle = %.3f (%.1f%%)", idx, dc, dc * 100)
        try:
            run_simulation_pipeline(
                duty_cycle=dc,
                ts=ts,
                qt_output_dir=qt_output_dir,
                experiment_id=idx
            )
        except Exception as e:
            LOGGER.exception("Experiment %d failed for D = %.4f: %s", idx, dc, e)
    LOGGER.info("PARAMETRIC RUNS ENDED")

# ----------------------------------------------------------------------
# Script Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
