#!/usr/bin/env python3
"""
Reaction Filtering Based on Thresholds

This script:
- Loads reaction rates over time.
- Computes the integral of each reaction.
- Filters them based on group-specific thresholds and pathway tags.
- Outputs both detailed and minimal filtered results.
"""

import numpy as np
import sys
import os
import logging
import json
from typing import Tuple, List
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
# Data Loaders
# ----------------------------------------------------------------------
def load_rates(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load time and rate data from Qt file.

    Returns:
        times (1D np.ndarray)
        rates (2D np.ndarray) shape = [timesteps x reactions]
    """
    try:
        data = np.loadtxt(file_path, skiprows=1)
        return data[:, 0], data[:, 1:]
    except Exception as e:
        logger.error(f"Failed to load rates file '{file_path}': {e}")
        raise

def load_reaction_names(file_path: str) -> List[str]:
    """
    Load reaction names from qt_reactions_list.txt (strip index, normalize).

    Returns:
        List of cleaned lowercase reaction strings.
    """
    try:
        with open(file_path, 'r') as f:
            return [line.strip().split(maxsplit=1)[1].replace(" ", "").lower() for line in f if line.strip()]
    except Exception as e:
        logger.error(f"Failed to load reaction names from '{file_path}': {e}")
        raise

# ----------------------------------------------------------------------
# Processing Functions
# ----------------------------------------------------------------------
def integrate_rates(times: np.ndarray, rates: np.ndarray) -> np.ndarray:
    """
    Compute area under the rate curve for each reaction.
    """
    return np.trapz(rates, times, axis=0)

def filter_by_group_thresholds(
    integrals: np.ndarray,
    reaction_names: List[str],
    reaction_groups: dict,
    group_thresholds: dict,
    pathway_tags: dict
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Apply threshold filtering by reaction group, preserving tagged reactions.

    Returns:
        low: indices of reactions below threshold
        high: indices of reactions to retain
        unprotected_discarded: names of discarded reactions not tagged
    """
    low, high = [], []
    unprotected_discarded = []

    reaction_to_group = {
        reaction.lower(): group
        for group, reactions in reaction_groups.items()
        for reaction in reactions
    }

    for idx, name in enumerate(reaction_names):
        group = reaction_to_group.get(name)
        threshold = group_thresholds.get(group)

        if group is None:
            logger.warning(f"Reaction not assigned to any group: {name}")
            continue
        if threshold is None:
            logger.warning(f"No threshold defined for group: {group}")
            continue

        if pathway_tags.get(name, False):
            high.append(idx)
        elif integrals[idx] < threshold:
            low.append(idx)
            if not pathway_tags.get(name, False):
                unprotected_discarded.append(name)
        else:
            high.append(idx)

    return np.array(low), np.array(high), unprotected_discarded

# ----------------------------------------------------------------------
# Output Writers
# ----------------------------------------------------------------------
def save_detailed_results(
    retained: np.ndarray,
    discarded: np.ndarray,
    integrals: np.ndarray,
    reaction_names: List[str],
    output_path: str = "reacciones_filtradas.txt"
):
    """
    Save human-readable summary of filtered reactions and their integrals.
    """
    try:
        with open(output_path, 'w') as f:
            f.write("# Retained reactions:\n")
            for idx in retained:
                f.write(f"R{idx+1:03d} - Integral: {integrals[idx]:.3e} - {reaction_names[idx]}\n")

            f.write("\n# Discarded reactions (below threshold):\n")
            for idx in discarded:
                f.write(f"R{idx+1:03d} - Integral: {integrals[idx]:.3e} - {reaction_names[idx]}\n")

        logger.info(f"Detailed reaction filter saved to {output_path}")
    except Exception as e:
        logger.error(f"Error writing detailed results to '{output_path}': {e}")
        raise

def save_retained_indices(
    retained: np.ndarray,
    reaction_names: List[str],
    output_path: str = "reacciones_usadas.txt"
):
    """
    Save the retained reaction indices in minimal format for downstream use.
    """
    try:
        with open(output_path, 'w') as f:
            for idx in retained:
                f.write(f"{idx+1} {reaction_names[idx]}\n")
        logger.info(f"Retained reaction indices saved to {output_path}")
    except Exception as e:
        logger.error(f"Error writing retained reactions to '{output_path}': {e}")
        raise

# ----------------------------------------------------------------------
# Main Filtering Function
# ----------------------------------------------------------------------
def run_filter_thresholds(active_thresholds_path: str) -> Tuple[bool, bool]:
    """
    Full filtering routine that:
    - Loads Qt rate and reaction data
    - Loads config files (groups, thresholds, tags)
    - Filters reactions
    - Saves outputs

    Returns:
        reactions_discarded (bool)
        all_protected (bool)
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    HELPER_FILES = PROJECT_ROOT / "src/helper_files"
    QT2_FILES = PROJECT_ROOT / "data/QtFiles/Qt_2"
    reaction_groups_path = HELPER_FILES / "reaction_groups.json"
    pathway_tags_path = HELPER_FILES / "reactions_tagged.json"
    rates_file = QT2_FILES / "qt_rates.txt"
    names_file = QT2_FILES / "qt_reactions_list.txt"


    try:
        logger.info("Loading data files...")
        times, rates = load_rates(rates_file)
        reaction_names = load_reaction_names(names_file)

        logger.info("Loading group thresholds...")
        with open(active_thresholds_path, 'r') as f:
            group_thresholds = json.load(f)

        logger.info("Loading reaction-to-group mapping...")
        with open(reaction_groups_path, 'r') as f:
            reaction_groups = json.load(f)

        logger.info("Loading reaction pathway tags...")
        with open(pathway_tags_path, 'r') as f:
            raw_tags = json.load(f)
            pathway_tags = {k.replace(" ", "").lower(): v for k, v in raw_tags.items()}

        logger.info("Computing integrals...")
        integrals = integrate_rates(times, rates)

        logger.info("Filtering by group-specific thresholds and pathway preservation...")
        low, high, unprotected_discarded = filter_by_group_thresholds(
            integrals, reaction_names, reaction_groups, group_thresholds, pathway_tags
        )

        logger.info(f"{len(low)} reactions discarded. {len(high)} retained.")
        save_detailed_results(high, low, integrals, reaction_names, HELPER_FILES / "reacciones_filtradas.txt")
        save_retained_indices(high, reaction_names, HELPER_FILES / "reacciones_usadas.txt")

        logger.info("Filtering complete.")
        reactions_discarded = len(low) > 0
        all_protected = reactions_discarded and len(unprotected_discarded) == 0

        return reactions_discarded, all_protected

    except Exception as e:
        logger.critical(f"Critical error during filtering: {e}")
        sys.exit(2)

# ----------------------------------------------------------------------
# CLI Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        logger.critical("Missing input: path to active_thresholds.json required.")
        sys.exit(1)

    active_thresholds_path = sys.argv[1]
    run_filter_thresholds(active_thresholds_path)
