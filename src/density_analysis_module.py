#!/usr/bin/env python3
"""
Density Analysis Module for ZDPlasKin Output

Computes:
  • Per-period ΔNO and Δ(N-radical) densities
  • Overall de-NOx efficiency

Outputs:
  • period_deltas.csv: Δ per species, per switching period
  • summary.txt: global efficiency and stats
  • PNG plots for each tracked species

Refactored as a reusable class: `DensityAnalyzer`
"""

import re
from pathlib import Path
import pandas as pd
import logging
from typing import Dict
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
# Utility Functions
# ----------------------------------------------------------------------
def parse_species_list(path: Path) -> Dict[str, str]:
    """
    Reads Qt species list and returns {species_id: name}.
    """
    id_to_name = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(None, 1)
            if len(parts) == 2:
                sid, name = parts
                id_to_name[sid] = name.strip()
    return id_to_name

def is_n_radical(name: str) -> bool:
    """
    Identify N* radical species (not NO or NO2).
    """
    return bool(re.match(r"^N(\(|\s|$|\^)", name)) and name not in {"NO", "NO2"}

# ----------------------------------------------------------------------
# Main Analyzer Class
# ----------------------------------------------------------------------
class DensityAnalyzer:
    """
    Analyzes ZDPlasKin output to compute de-NOx efficiency and per-period deltas.

    Attributes:
        qt_dir: directory containing qt_species_list.txt and qt_densities.txt
        ts: switching period in seconds
        output_dir: directory to write plots and output files
        sim_time: optional simulation runtime (for summary.txt)
    """

    def __init__(self, qt_dir: Path, ts: float, output_dir: Path, sim_time: float):
        self.qt_dir = qt_dir
        self.ts = ts
        self.output_dir = output_dir
        self.sim_time = sim_time
        self.sp_file = qt_dir / "qt_species_list.txt"
        self.dens_file = qt_dir / "qt_densities.txt"

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_deltas(self, merged_df: pd.DataFrame, name: str):
        """
        Save time-series plot of delta_{name}.
        """
        x = merged_df["end_time"]
        plt.figure(figsize=(8, 5))
        plt.plot(x, merged_df[f"delta_{name}"], label=f"Δ{name}")
        plt.title(f"Variation of {name}", fontsize=14)
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Δ Density (cm⁻³)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(self.output_dir / f"delta_{name}.png", dpi=300)
        plt.close()
        LOGGER.info(f"Saved plot for Δ{name}")

    def run(self) -> float:
        """
        Executes full density analysis and saves:
          - period_deltas.csv
          - summary.txt
          - one plot per species

        Returns:
            float: de-NOx efficiency.
        """
        if not self.sp_file.exists() or not self.dens_file.exists():
            raise FileNotFoundError("Missing qt_species_list.txt or qt_densities.txt")

        id_to_name = parse_species_list(self.sp_file)

        # Identify species
        id_NO = next((sid for sid, n in id_to_name.items() if n.strip() == "NO"), None)
        id_NO2 = next((sid for sid, n in id_to_name.items() if n.strip() == "NO2"), None)
        if id_NO is None or id_NO2 is None:
            raise ValueError("NO or NO2 not found in species list")

        n_radical_ids = [sid for sid, name in id_to_name.items() if is_n_radical(name)]
        if not n_radical_ids:
            raise ValueError("No N-radical species found")

        # Load density file
        df = pd.read_csv(self.dens_file, delim_whitespace=True, comment="#")
        df.columns = [str(c) for c in df.columns]

        time_col = "Time_s"
        if time_col not in df.columns:
            raise ValueError(f"Expected column '{time_col}' in qt_densities")

        # Partition by switching periods
        df["period"] = (df[time_col] // self.ts).astype(int)
        first_idx = df.groupby("period").head(1).index
        last_idx = df.groupby("period").tail(1).index

        # Extract beginning and end densities
        start_df = df.loc[first_idx, ["period", time_col, id_NO] + n_radical_ids]
        end_df   = df.loc[last_idx,  ["period", time_col, id_NO] + n_radical_ids]

        # Rename for merging
        start_df = start_df.rename(columns={
            time_col: "start_time",
            id_NO: "NO_start",
            **{sid: f"{id_to_name[sid]}_start" for sid in n_radical_ids}
        })
        end_df = end_df.rename(columns={
            time_col: "end_time",
            id_NO: "NO_end",
            **{sid: f"{id_to_name[sid]}_end" for sid in n_radical_ids}
        })

        # Merge and compute deltas
        merged = pd.merge(start_df, end_df, on="period")
        merged["delta_NO"] = merged["NO_end"] - merged["NO_start"]
        for sid in n_radical_ids:
            sname = id_to_name[sid]
            merged[f"delta_{sname}"] = merged[f"{sname}_end"] - merged[f"{sname}_start"]

        delta_cols = ["delta_NO"] + [f"delta_{id_to_name[sid]}" for sid in n_radical_ids]
        out_cols = ["period", "start_time", "end_time"] + [f"{c}_scaled" for c in delta_cols]

        deltas_path = self.output_dir / "period_deltas.csv"
        summary_path = self.output_dir / "summary.txt"

        SCALING_FACTOR = 1  # If needed to rescale densities

        for col in delta_cols:
            merged[f"{col}_scaled"] = merged[col] * SCALING_FACTOR

        merged[out_cols].to_csv(deltas_path, index=False, float_format="%.5e")
        LOGGER.info("Saved per-period deltas to %s", deltas_path)

        # Efficiency = 1 - (NO + NO2 final) / NO initial
        NO_initial = df.loc[0, id_NO]
        NO_final = df.loc[df.index[-1], id_NO]
        NO2_final = df.loc[df.index[-1], id_NO2]
        eff = 1.0 - (NO_final + NO2_final) / NO_initial

        # Summary report
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(f"de-NOx efficiency      : {eff:.6f}\n")
            f.write(f"Initial NO density     : {float(NO_initial):.5e}  [cm⁻³]\n")
            f.write(f"Final   NO density     : {float(NO_final):.5e}\n")
            f.write(f"Final  NO2 density     : {float(NO2_final):.5e}\n")
            f.write(f"Switch-period Ts       : {self.ts:.3e} s\n")
            f.write(f"Simulation time (s)    : {self.sim_time:.3f}\n" if self.sim_time else "Simulation time (s)    : Not available\n")

        for sid in delta_cols:
            self.plot_deltas(merged, sid)

        LOGGER.info("Computed de-NOx efficiency: %.6f", eff)
        return eff
