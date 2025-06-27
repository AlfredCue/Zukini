#!/usr/bin/env python3
"""
Sensitivity-Based Reaction Protection

This script computes a stoichiometric sensitivity matrix from QtPlasKin output,
identifying reactions that significantly contribute to target species dynamics.

Those reactions are marked as protected in `reactions_tagged.json`.
"""

import numpy as np
import pandas as pd
import json
import sys
import os
from typing import List, Tuple
from pathlib import Path

# ----------------------------------------------------------------------
# Data Loading Utilities
# ----------------------------------------------------------------------

def load_species_list(filepath: str) -> pd.DataFrame:
    """
    Load QtPlasKin species index + names from txt file.

    Returns:
        DataFrame with columns: Index, Species
    """
    return pd.read_csv(filepath, sep=r'\s+', header=None, names=["Index", "Species"])

def load_reactions_list(filepath: str) -> pd.DataFrame:
    """
    Load QtPlasKin reaction index + formulas from txt file.

    Returns:
        DataFrame with columns: Index, Reaction
    """
    return pd.read_csv(filepath, sep=r'\s+', header=None, names=["Index", "Reaction"], engine='python')

def load_stoichiometry_matrix(filepath: str) -> np.ndarray:
    """
    Load stoichiometry matrix from qt_matrix.txt.

    Returns:
        Matrix of shape [species × reactions]
    """
    return np.loadtxt(filepath)

def load_rates(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load reaction rate time series from qt_rates.txt.

    Returns:
        time (np.ndarray), rates (reactions × time_steps)
    """
    df = pd.read_csv(filepath, sep=r'\s+')
    time = df.iloc[:, 0].values
    rates = df.iloc[:, 1:].values
    return time, rates.T  # transpose so shape = [reactions x time]

# ----------------------------------------------------------------------
# Core Computational Steps
# ----------------------------------------------------------------------

def compute_species_source_terms(stoich_matrix: np.ndarray, rates_array: np.ndarray) -> np.ndarray:
    """
    Compute Qij(t): contribution of reaction j to species i at each time step.

    Returns:
        3D array [species, reactions, time_steps]
    """
    species_count, reaction_count = stoich_matrix.shape
    reaction_count2, time_steps = rates_array.shape
    assert reaction_count == reaction_count2, "Mismatch between stoich and rate data."

    Qij_time = np.zeros((species_count, reaction_count, time_steps))
    for j in range(reaction_count):
        Qij_time[:, j, :] = np.outer(stoich_matrix[:, j], rates_array[j, :])

    return Qij_time

def integrate_over_time(time: np.ndarray, Qij_time: np.ndarray) -> np.ndarray:
    """
    Integrate absolute value of Qij(t) over time.

    Returns:
        2D array [species x reactions]
    """
    abs_Qij_time = np.abs(Qij_time)
    dt = np.diff(time, prepend=time[0])
    integral = np.sum(abs_Qij_time * dt, axis=2)
    return integral

def compute_sensitivity_matrix(integrated_Qij: np.ndarray, integrated_Qi: np.ndarray) -> np.ndarray:
    """
    Compute normalized sensitivity: S[i,j] = Qij / Qi

    Avoids divide-by-zero errors.

    Returns:
        Sensitivity matrix [species x reactions]
    """
    integrated_Qi[integrated_Qi == 0] = 1e-30
    return integrated_Qij / integrated_Qi[:, np.newaxis]

# ----------------------------------------------------------------------
# Main Analysis & Update
# ----------------------------------------------------------------------

def run_sensitivity_filter(target_species: List[str]) -> None:
    """
    Executes full sensitivity analysis pipeline and updates reactions_tagged.json
    to protect reactions affecting target_species above a threshold.
    """
    # ==== CONFIG ====
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    QT1_FILES = PROJECT_ROOT / "Data/QtFiles/Qt_1"
    species_list_file = QT1_FILES / "qt_species_list.txt"
    reactions_list_file = QT1_FILES / "qt_reactions_list.txt"
    stoich_matrix_file = QT1_FILES / "qt_matrix.txt"
    rates_file = QT1_FILES / "qt_rates.txt"
    reactions_tagged_file = PROJECT_ROOT / "src/helper_files/reactions_tagged.json"
    sensitivity_threshold = 1e-30

    # ==== LOAD ====
    print("[INFO] Loading QtPlasKin outputs...")
    species_df = load_species_list(species_list_file)
    reactions_df = load_reactions_list(reactions_list_file)
    stoich_matrix = load_stoichiometry_matrix(stoich_matrix_file)
    time_array, rates_array = load_rates(rates_file)

    # ==== CALCULATIONS ====
    print("[INFO] Computing Qij contributions...")
    Qij_time = compute_species_source_terms(stoich_matrix, rates_array)

    print("[INFO] Integrating Qij over time...")
    integrated_Qij = integrate_over_time(time_array, Qij_time)
    integrated_Qi = np.sum(np.abs(integrated_Qij), axis=1)

    print("[INFO] Computing sensitivity matrix...")
    sensitivity_matrix = compute_sensitivity_matrix(integrated_Qij, integrated_Qi)

    # ==== DETERMINE REACTIONS TO PROTECT ====
    species_to_index = {row["Species"]: row["Index"] - 1 for _, row in species_df.iterrows()}
    protected_reactions = set()

    for specie in target_species:
        if specie not in species_to_index:
            print(f"[WARNING] Species '{specie}' not found in species list!")
            continue
        idx = species_to_index[specie]
        sensitive_idxs = np.where(sensitivity_matrix[idx, :] >= sensitivity_threshold)[0]
        protected_reactions.update(sensitive_idxs)

    print(f"[INFO] {len(protected_reactions)} reactions will be protected based on sensitivity.")

    # ==== UPDATE TAGGED FILE ====
    print("[INFO] Updating reactions_tagged.json...")
    with open(reactions_tagged_file, 'r') as f:
        tagged_data = json.load(f)

    reaction_names = reactions_df["Reaction"].values
    updated_count = 0

    for reaction_idx in protected_reactions:
        if 0 <= reaction_idx < len(reaction_names):
            formula = reaction_names[reaction_idx]
            if formula in tagged_data:
                tagged_data[formula] = True
                updated_count += 1
            else:
                print(f"[WARNING] Formula not found in reactions_tagged.json: {formula}")

    with open(reactions_tagged_file, 'w') as f:
        json.dump(tagged_data, f, indent=4)

    print(f"[INFO] Successfully updated {updated_count} reactions in reactions_tagged.json.")

# ----------------------------------------------------------------------
# CLI Entry Point
# ----------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("[ERROR] No species provided as input argument.")
        print("Usage: python sensitivity_filter.py NO,NO2,N")
        sys.exit(1)

    species_input = sys.argv[1]
    target_species = species_input.split(",")

    run_sensitivity_filter(target_species)
