#!/usr/bin/env python3
"""
Brute Force Recovery Module

This module decides whether to reintroduce previously filtered reactions
based on increased simulation time (runtime regression). If the simulation
is slower after filtering, reactions are marked back as protected.

Also updates:
- reactions_tagged.json
- recovered_reactions.json
"""

import json
import re
import logging
from pathlib import Path
from typing import Optional, Union
from datetime import datetime

# ----------------------------------------------------------------------
# Regex to Extract Simulation Time from ZDPlasKin Logs
# ----------------------------------------------------------------------
_SIM_TIME_RE = re.compile(
    r"simulation\s*time\s*:?\s*([0-9.+\-eE]+)",
    re.I  # case-insensitive
)

_SENS_RE = re.compile(
    r"NO2\s*selectivity\s*:?\s*([0-9.+\-eE]+)",
    re.I
)

#-----------------------------------------------------------------------
# Logging setup
#-----------------------------------------------------------------------
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
# Simulation Time Extraction
# ----------------------------------------------------------------------
def extract_simulation_time(log_path: Union[str,Path]) -> Optional[float]:
    """
    Parses simulation time from a ZDPlasKin output log file.

    Args:
        log_path: Path to log file.

    Returns:
        Float (seconds) if found, None otherwise.
    """
    log_path = Path(log_path)
    if not log_path.exists():
        LOGGER.error(f"File not found: {log_path}")
        return None

    with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            match = _SIM_TIME_RE.search(line)
            if match:
                try:
                    return float(match.group(1))
                except ValueError as exc:
                    LOGGER.error(f"Could not parse '{match.group(1)}' as float: {exc}")
                    return None

    LOGGER.error("Simulation time not found in log.")
    return None

def extract_selectivity(log_path: Union[str,Path])-> Optional[float]:
    log_path = Path(log_path)
    if not log_path.exists():
        LOGGER.error(f"File not found: {log_path}")
        return None
    with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            match = _SENS_RE.search(line)
            if match:
                try:
                    return float(match.group(1))
                except ValueError as exc:
                    LOGGER.error(f"Could not parse '{match.group(1)}' as float: {exc}")
                    return None

    LOGGER.error("NO2 Selectivity not found in log.")
    return None


# ----------------------------------------------------------------------
# Parse Discarded Reactions
# ----------------------------------------------------------------------
def parse_filtered_reactions(filtered_filepath):
    """
    Parses the reaction numbers marked as discarded in 'reacciones_filtradas.txt'.

    Returns:
        List[int] of reaction numbers.
    """
    with open(filtered_filepath, 'r') as f:
        lines = f.readlines()

    discarded = []
    reading_discarded = False

    for line in lines:
        line = line.strip()
        if line.startswith("# Discarded reactions"):
            reading_discarded = True
            continue
        if reading_discarded and line:
            match = re.match(r"R(\d+)", line)
            if match:
                discarded.append(int(match.group(1)))

    return discarded

# ----------------------------------------------------------------------
# Load Reaction Index → Formula Mapping
# ----------------------------------------------------------------------
def load_reactions_list(filepath):
    """
    Loads reaction index-to-formula mapping from qt_reactions_list.txt.

    Returns:
        Dict[int, str]
    """
    reaction_map = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                idx, formula = parts
                reaction_map[int(idx)] = formula
    return reaction_map

# ----------------------------------------------------------------------
# Main Brute Force Recovery Routine
# ----------------------------------------------------------------------
def brute_force_recovery(
    reference_time_sec,
    reduced_log,
    filtered_reactions_file,
    reactions_tagged_file,
    force,
    recovered_output_file="recovered_reactions.json"
) -> bool:
    """
    Decide whether to recover reactions based on increased runtime.

    Args:
        reference_time_sec (float): Time of baseline simulation (original model).
        reduced_log (str): Path to reduced simulation log.
        filtered_reactions_file (str): File with discarded reaction IDs.
        reactions_tagged_file (str): File to mark reactions as "protected" (True).
        force (bool): Force recovery regardless of time check.
        recovered_output_file (str): File to save list of recovered reaction formulas.

    Returns:
        True if reactions were recovered, False otherwise.
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    QT1_FILES = PROJECT_ROOT / "Data/QtFiles/Qt_1"
    
    LOGGER.info("Extracting reduced simulation time...")
    reduced_time = extract_simulation_time(reduced_log)

    if reduced_time is None:
        LOGGER.error("Could not extract reduced simulation time.")
        return False

    LOGGER.info(f"Reference simulation time: {reference_time_sec:.4f} s")
    LOGGER.info(f"Reduced simulation time: {reduced_time:.4f} s")

    if reduced_time > reference_time_sec * 1.1 or force:
        LOGGER.warning("Reduced simulation is slower or forced. Recovering reactions...")
        recovered_formulas = []

        # Step 1: Parse discarded reaction indices
        discarded_reactions = parse_filtered_reactions(filtered_reactions_file)

        if not discarded_reactions:
            LOGGER.info("No discarded reactions found. Nothing to recover.")
            return False

        LOGGER.info(f"{len(discarded_reactions)} discarded reactions will be recovered.")

        # Step 2: Load reaction index → formula
        reaction_map = load_reactions_list(QT1_FILES / "qt_reactions_list.txt")

        # Step 3: Load current reaction tags
        with open(reactions_tagged_file, 'r') as f:
            tagged_data = json.load(f)

        # Step 4: Update tags
        recovered_count = 0
        for reaction_number in discarded_reactions:
            formula = reaction_map.get(reaction_number)
            if formula and formula in tagged_data:
                tagged_data[formula] = True
                recovered_formulas.append(formula)
                recovered_count += 1
            else:
                LOGGER.warning(f"Reaction R{reaction_number:03d} not found in reactions_tagged.json")

        # Step 5: Save updated tags
        with open(reactions_tagged_file, 'w') as f:
            json.dump(tagged_data, f, indent=4)

        LOGGER.info(f"Updated {recovered_count} reactions in {reactions_tagged_file}")

        # Step 6: Save recovered list
        with open(recovered_output_file, 'w') as f:
            json.dump(recovered_formulas, f, indent=4)

        LOGGER.info(f"Recovered reaction formulas saved to {recovered_output_file}")
        return True

    else:
        LOGGER.info("Reduced simulation is faster. No recovery needed.")
        return False

# ----------------------------------------------------------------------
# CLI Manual Test Mode
# ----------------------------------------------------------------------
if __name__ == "__main__":
    recovered = brute_force_recovery(
        reference_time_sec=692.0,
        reduced_log="ZDPlaskin/simulation_output.txt",
        filtered_reactions_file="reacciones_filtradas.txt",
        reactions_tagged_file="reactions_tagged.json",
        force=False
    )

    if recovered:
        print("[FINAL INFO] Some reactions were recovered based on simulation time check.")
    else:
        print("[FINAL INFO] No reactions were recovered. Reduced model is stable.")
