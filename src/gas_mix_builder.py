#!/usr/bin/env python3
"""
Gas Mixture Builder

Updates NO and O2 ppm values in a Fortran source file (e.g. dbd0D.F90)
based on concentrations provided via a JSON configuration.

Intended for automated modification of simulation conditions in parametric studies.
"""

import json
import re
import logging
from pathlib import Path
from typing import Union

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
# GasMixture Updater Class
# ----------------------------------------------------------------------
class GasMixBuilder:
    """
    Replaces NO and O2 concentration constants in a Fortran file using values
    provided in a JSON file.

    Attributes:
        json_file: Path to JSON with {"NO": value, "O2": value}
        fortran_input: Path to the source .F90 file to read from
        fortran_output: Path to write the modified .F90 file
    """

    def __init__(self, json_file: Union[str, Path], fortran_input: Union[str, Path], fortran_output: Union[str, Path]):
        self.json_file = Path(json_file)
        self.fortran_input = Path(fortran_input)
        self.fortran_output = Path(fortran_output)

    def update(self) -> None:
        """
        Performs the actual substitution of gas concentrations in the Fortran file.
        """
        if not self.json_file.exists():
            raise FileNotFoundError(f"JSON file not found: {self.json_file}")
        if not self.fortran_input.exists():
            raise FileNotFoundError(f"Fortran input file not found: {self.fortran_input}")

        # Load gas concentrations
        with open(self.json_file, "r", encoding="utf-8") as f:
            gas = json.load(f)

        try:
            no_value = float(gas["NO"])
            o2_value = float(gas["O2"])
        except KeyError as err:
            raise KeyError("JSON must contain both 'NO' and 'O2' keys.") from err

        # Define matching patterns for Fortran parameters
        num_pattern = r"[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?"
        no_pattern = re.compile(rf"(real\(8\),\s*parameter\s*::\s*no_ppm\s*=\s*)({num_pattern})", re.I)
        o2_pattern = re.compile(rf"(real\(8\),\s*parameter\s*::\s*o2_ppm\s*=\s*)({num_pattern})", re.I)

        # Read and modify Fortran source
        with open(self.fortran_input, "r", encoding="utf-8") as f:
            source = f.read()

        new_source = no_pattern.sub(rf"\g<1>{no_value:g}", source)
        new_source = o2_pattern.sub(rf"\g<1>{o2_value:g}", new_source)

        # Write to output (can be same file)
        with open(self.fortran_output, "w", encoding="utf-8") as f:
            f.write(new_source)

        LOGGER.info("Gas concentrations updated and written to: %s", self.fortran_output.resolve())

# ----------------------------------------------------------------------
# CLI/Dev Test Mode
# ----------------------------------------------------------------------
if __name__ == "__main__":
    builder = GasMixBuilder(
        json_file="gas_mixture.json",
        fortran_input="ZDPlaskin/dbd0D.F90",
        fortran_output="ZDPlaskin/dbd0D.F90"
    )
    builder.update()
