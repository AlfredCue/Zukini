#!/usr/bin/env python3
"""
Robustness-Analysis Results Aggregator

This script scans a directory tree containing robustness analysis results
and compiles per-run metrics into a single CSV file.

Expected folder structure:
    robustness_analysis_results/
      └─ T300_G_NO_500_O2_20_Dc_0_1/
         └─ filtered1/
            └─ metrics_summary.txt
         └─ filtered2/
            └─ metrics_summary.txt

Output:
    A CSV with columns like temperature, NO/O2 ppm, duty, filtered set ID,
    and metrics (efficiency, trust score, runtimes).

Usage:
    python run_robustness_check.py --root robustness_analysis_results --out compiled_results.csv
"""

import re
import argparse
import logging
from pathlib import Path
from typing import Optional, Dict, List, Union

import pandas as pd

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
LOGGER = logging.getLogger("results-aggregator")

# ----------------------------------------------------------------------
# Regex Patterns
# ----------------------------------------------------------------------
RUN_RE = re.compile(
    r"""
    ^T(?P<temp>\d+)              _        # temperature
    G_NO_(?P<no>\d+)             _        # NO ppm
    O2_(?P<o2>\d+)               _        # O2 ppm
    Dc_(?P<dc>[0-9_]+)$                   # duty cycle (with underscores)
    """, re.X
)

NUM_RE = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")


# ----------------------------------------------------------------------
# Aggregator Class
# ----------------------------------------------------------------------
class RobustnessResultsAggregator:
    """
    Scans a directory tree for `metrics_summary.txt` files and compiles them
    into a single CSV DataFrame.
    """

    def __init__(self):
        self.run_re = RUN_RE
        self.num_re = NUM_RE

    def compile(
        self,
        root: Union[str, Path] = "robustness_analysis_results",
        out_csv: Union[str, Path] = "compiled_results.csv",
    ) -> pd.DataFrame:
        """
        Main method: scans the directory tree and generates a DataFrame.

        Args:
            root: Root directory to scan.
            out_csv: Output CSV file.

        Returns:
            Compiled DataFrame.
        """
        root = Path(root)
        out_csv = Path(out_csv)

        LOGGER.info("Scanning for metrics_summary.txt under %s …", root.resolve())
        records: List[Dict] = []

        for metrics_path in root.rglob("metrics_summary.txt"):
            try:
                filtered_dir = metrics_path.parent
                run_dir = filtered_dir.parent
                run_name = run_dir.name

                match = self.run_re.match(run_name)
                if not match:
                    LOGGER.warning("Skipping unrecognized folder name: %s", run_name)
                    continue

                # Parse metadata from folder name
                temp = int(match.group("temp"))
                NO_ppm = int(match.group("no"))
                O2_ppm = int(match.group("o2"))
                duty = float(match.group("dc").replace("_", "."))

                tag = filtered_dir.name.lower()
                if not tag.startswith("filtered"):
                    LOGGER.warning("Skipping path without filteredX folder: %s", metrics_path)
                    continue

                filtered_set = int(tag.replace("filtered", "").strip())
                metrics = self._parse_metrics_file(metrics_path)

                records.append({
                    "temp": temp,
                    "NO_ppm": NO_ppm,
                    "O2_ppm": O2_ppm,
                    "duty_cycle": duty,
                    "filtered_set": filtered_set,
                    **metrics,
                })
            except Exception as exc:
                LOGGER.exception("Failed to parse %s: %s", metrics_path, exc)

        if not records:
            raise RuntimeError("No metrics files found — nothing to compile.")

        df = (
            pd.DataFrame.from_records(records)
            .sort_values(["temp", "NO_ppm", "O2_ppm", "duty_cycle", "filtered_set"])
            .reset_index(drop=True)
        )

        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        LOGGER.info("Compiled %d rows → %s", len(df), out_csv.resolve())
        return df

    # ------------------------------------------------------------------
    # Internal Helpers
    # ------------------------------------------------------------------
    def _extract_number(self, line: str) -> Optional[float]:
        match = self.num_re.search(line)
        return float(match.group(1)) if match else None

    def _parse_metrics_file(self, path: Path) -> Dict[str, Optional[float]]:
        vals = {
            "eff_full": None,
            "eff_reduced": None,
            "delta_eff": None,
            "trust_score": None,
            "sim_time_full": None,
            "sim_time_red": None,
        }

        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                if "de-NOx Efficiency (full)" in line:
                    vals["eff_full"] = self._extract_number(line)
                elif "de-NOx Efficiency (reduced)" in line:
                    vals["eff_reduced"] = self._extract_number(line)
                elif "de-NOx Δ Efficiency" in line:
                    vals["delta_eff"] = self._extract_number(line)
                elif "Overall Trust Score" in line:
                    vals["trust_score"] = self._extract_number(line)
                elif "Simulation Time (s) full" in line:
                    vals["sim_time_full"] = self._extract_number(line)
                elif "Simulation Time (s) reduced" in line:
                    vals["sim_time_red"] = self._extract_number(line)

        return vals


# ----------------------------------------------------------------------
# CLI Interface
# ----------------------------------------------------------------------
def _cli():
    parser = argparse.ArgumentParser(description="Compile robustness-analysis metrics.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("robustness_analysis_results"),
        help="Root directory containing all simulation folders."
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("compiled_results.csv"),
        help="Output CSV path."
    )
    args = parser.parse_args()

    aggregator = RobustnessResultsAggregator()
    aggregator.compile(args.root, args.out)


if __name__ == "__main__":
    _cli()
