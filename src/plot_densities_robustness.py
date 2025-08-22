#!/usr/bin/env python3
"""
ZDPlasKin Robustness Analyzer

Compares full vs. reduced simulation outputs:
- Computes per-species trust metrics (R² + NRMSE-based)
- Generates plots for each species
- Aggregates a global weighted trust score
- Optionally compares de-NOx efficiencies
"""

from __future__ import annotations
import logging
from pathlib import Path
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt



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
# Data Loading Utilities
# ----------------------------------------------------------------------
def load_densities(file_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path, skiprows=1)
    return data[:, 0], data[:, 1:]

def load_species_names(file_path: Path) -> List[str]:
    with open(file_path, "r", encoding="utf-8") as f:
        return [line.strip().split(maxsplit=1)[1] for line in f if line.strip()]

def extract_density(species: str, species_names: List[str], density_matrix: np.ndarray) -> np.ndarray:
    try:
        idx = species_names.index(species)
    except ValueError as exc:
        raise ValueError(f"Species '{species}' not found in species list.") from exc
    return density_matrix[:, idx]

def extract_last_value(species: str, species_names: List[str], density_matrix: np.ndarray) -> np.ndarray:
    try:
        idx = species_names.index(species)
    except ValueError as exc:
        raise ValueError(f"Species '{species}' not found in species list.") from exc
    return density_matrix[-1, idx]

# ----------------------------------------------------------------------
# Metrics Calculations
# ----------------------------------------------------------------------
def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 0.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot

def nrmse_peak(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    rmse = np.sqrt(np.mean((y_true - y_pred) ** 2))
    peak = np.max(np.abs(y_true))
    return rmse / peak if peak > 0 else rmse

def trust_score(r2: float, nrmse: float, w_r2: float = 0.6, beta: float = 2.0) -> float:
    penalty = float(np.exp(-beta * nrmse))
    return np.clip(w_r2 * r2 + (1 - w_r2) * penalty, 0.0, 1.0)

# ----------------------------------------------------------------------
# Plotting Utilities
# ----------------------------------------------------------------------
def _style():
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 9,
        "ytick.labelsize": 11,
        "legend.fontsize": 12,
        "font.family": "sans-serif",
    })

def plot_comparison(species: str, time: np.ndarray, full: np.ndarray, reduced: np.ndarray,
                    trust_pct: float, out_dir: Path):
    _style()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time, full, label="Full", linewidth=2)
    ax.plot(time, reduced, label="Reduced", linewidth=2)
    ax.set_title(f"Density Evolution {species}", fontweight="bold")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Density (cm⁻³)")
    ax.text(0.02, 0.95, f"Trust: {trust_pct:.1f}%", transform=ax.transAxes,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"comparison_{species}.png")
    plt.close(fig)

def plot_bar(labels: List[str], values: np.ndarray, title: str, ylabel: str, colour: str, out_path: Path):
    _style()
    fig, ax = plt.subplots(figsize=(max(10, 0.35 * len(labels)), 6))
    bars = ax.bar(labels, values, color=colour, edgecolor="black")
    ax.set_title(title, fontweight="bold")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, values.max() * 1.15)
    ax.yaxis.grid(True, linestyle="--", alpha=0.6)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 0.5, f"{height:.1f}",
                ha="center", va="bottom", fontsize=8, rotation=60)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_path, dpi=600)
    plt.close(fig)

# ----------------------------------------------------------------------
# Analyzer Class
# ----------------------------------------------------------------------
class ZDPlasKinAnalyzer:
    def __init__(self, w_r2=0.6, beta=2.0, priority_boost=10.0):
        self.w_r2 = w_r2
        self.beta = beta
        self.priority_boost = priority_boost
        PROJECT_ROOT = Path(__file__).resolve().parents[1]
        self.IMAGES_DIR = PROJECT_ROOT / "images"
        self.OUTPUT_DIR = PROJECT_ROOT / "output/robustness_analysis_results"

    def analyze(self, qt_full: Path, qt_reduced: Path,
                priority_species: List[str],
                sim_time_reduced: float,
                sim_time_full: float,
                selectivity_full: float,
                selectivity_reduced: float,
                o3_full: float,
                o3_reduced: float) -> float:
        """
        Compare full and reduced Qt outputs and compute global trust metric.

        Outputs:
            - PNG plots per species
            - metrics_summary.txt
        Returns:
            float: weighted global trust score [0-100]
        """

        LOGGER.info("Loading Qt data …")
        t_full, dens_full = load_densities(qt_full / "qt_densities.txt")
        species_full = load_species_names(qt_full / "qt_species_list.txt")

        t_red, dens_red = load_densities(qt_reduced / "qt_densities.txt")
        species_red = load_species_names(qt_reduced / "qt_species_list.txt")

        if species_full != species_red:
            raise ValueError("Species ordering mismatch between full and reduced runs.")

        # Interpolate if time bases differ
        if not np.allclose(t_full, t_red, rtol=1e-3, atol=1e-3):
            dens_red_interp = np.array([
                np.interp(t_full, t_red, dens_red[:, i]) for i in range(dens_red.shape[1])
            ]).T
        else:
            dens_red_interp = dens_red

        trust_scores = []
        for i, sp in enumerate(species_full):
            y_true = dens_full[:, i]
            y_pred = dens_red_interp[:, i]
            r2 = r2_score(y_true, y_pred)
            nrmse = nrmse_peak(y_true, y_pred)
            trust = trust_score(r2, nrmse, w_r2=self.w_r2, beta=self.beta)
            trust_scores.append(trust * 100)
            plot_comparison(sp, t_full, y_true, y_pred, trust * 100, out_dir=self.IMAGES_DIR / "comparisons")

        # Weighted average trust
        trust_scores_arr = np.array(trust_scores)
        dt = np.mean(np.diff(t_full))
        participation = np.sum(dens_full, axis=0) * dt
        participation_pct = participation / participation.sum() * 100

        weights = []
        for p, name in zip(participation_pct, species_full):
            weight = p / 100.0
            if name.upper() in [s.upper() for s in priority_species]:
                weight *= self.priority_boost
            weights.append(weight)

        overall = float(np.sum(np.array(weights) * trust_scores_arr) / sum(weights))
        LOGGER.info("Overall weighted trust: %.2f", overall)

        # Efficiency metrics (if NO/NO2 available)
        try:
            id_no = species_full.index("NO")
            id_no2 = species_full.index("NO2")
            eff_full = 1.0 - (dens_full[-1, id_no] + dens_full[-1, id_no2]) / dens_full[0, id_no]
            eff_red  = 1.0 - (dens_red_interp[-1, id_no] + dens_red_interp[-1, id_no2]) / dens_red_interp[0, id_no]
            delta_eff = eff_red - eff_full
        except ValueError:
            LOGGER.warning("NO/NO2 not found. Efficiency skipped.")
            eff_full = eff_red = delta_eff = None

        # Save summary
        summary_path = self.OUTPUT_DIR / "metrics_summary.txt"
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(f"de-NOx Efficiency (full):     {eff_full:.4f}\n" if eff_full is not None else "N/A\n")
            f.write(f"de-NOx Efficiency (reduced):  {eff_red:.4f}\n" if eff_red is not None else "N/A\n")
            f.write(f"de-NOx Δ Efficiency:          {delta_eff:+.4f}\n\n" if delta_eff is not None else "N/A\n")
            f.write(f"Overall Trust Score:          {overall:.2f} %\n\n")
            f.write(f"Simulation Time (s) full:     {sim_time_full:.3f}\n")
            f.write(f"Simulation Time (s) reduced:  {sim_time_reduced:.3f}\n")
            f.write(f"Selectivity full:             {selectivity_full:.4f}\n")
            f.write(f"Selectivity reduced:          {selectivity_reduced:.4f}\n")
            f.write(f"Final O3 full:                {o3_full:.10f}\n")
            f.write(f"Final O3 reduced:             {o3_reduced:.10f}\n")
        return overall

    def plot_single_species(self, species: str, qt_full: Path, qt_reduced: Path) -> None:
        """
        Plot full vs reduced densities for a single species only.

        Saves:
            PNG under output_dir / "comparisons"
        """
        t_full, dens_full = load_densities(qt_full / "qt_densities.txt")
        t_red, dens_red = load_densities(qt_reduced / "qt_densities.txt")
        species_names = load_species_names(qt_full / "qt_species_list.txt")

        # Interpolate reduced if necessary
        if not np.allclose(t_full, t_red, rtol=1e-3, atol=1e-3):
            dens_red = np.vstack([
                np.interp(t_full, t_red, dens_red[:, i])
                for i in range(dens_red.shape[1])
            ]).T

        y_true = extract_density(species, species_names, dens_full)
        y_pred = extract_density(species, species_names, dens_red)

        trust_pct = trust_score(r2_score(y_true, y_pred), nrmse_peak(y_true, y_pred), self.w_r2, self.beta) * 100.0

        plot_comparison(species, t_full, y_true, y_pred, trust_pct,out_dir= self.IMAGES_DIR / "comparisons")
