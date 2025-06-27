#!/usr/bin/env python3
"""
Constant Temperature Builder

Updates the `gas_temperature` value in a Fortran input file (dbd0D.F90)
based on a JSON configuration. Modularized for automated parametric sweeps.
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
# Temperature Updater Class
# ----------------------------------------------------------------------
class TemperatureBuilder:
    """
    Injects a new gas temperature value into a Fortran source file,
    replacing the constant definition.

    Attributes:
        json_file: Path to config with {"gas_temperature": <value>}
        fortran_input: Source .F90 file to read from
        fortran_output: Path to write modified version (can be same)
    """

    def __init__(self, json_file: Union[str, Path], fortran_input: Union[str, Path], fortran_output: Union[str, Path]):
        self.json_file = Path(json_file)
        self.fortran_input = Path(fortran_input)
        self.fortran_output = Path(fortran_output)

    def update(self) -> None:
        """
        Parses the JSON and replaces the `gas_temperature` constant in the Fortran source.
        """
        if not self.json_file.exists():
            raise FileNotFoundError(f"JSON file not found: {self.json_file}")
        if not self.fortran_input.exists():
            raise FileNotFoundError(f"Fortran input file not found: {self.fortran_input}")

        # Load temperature value
        with open(self.json_file, "r", encoding="utf-8") as f:
            config = json.load(f)

        try:
            temperature = float(config["gas_temperature"])
        except KeyError as err:
            raise KeyError("JSON must contain the 'gas_temperature' key.") from err

        # Regex pattern to match constant assignment
        num_pattern = r"[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?"
        temp_pattern = re.compile(rf"(real\(8\),\s*parameter\s*::\s*gas_temperature\s*=\s*)({num_pattern})", re.I)

        # Modify source file
        with open(self.fortran_input, "r", encoding="utf-8") as f:
            source = f.read()

        new_source = temp_pattern.sub(rf"\g<1>{temperature:g}", source)

        # Write to output
        with open(self.fortran_output, "w", encoding="utf-8") as f:
            f.write(new_source)

        LOGGER.info("Gas temperature updated to %.2f K → %s", temperature, self.fortran_output.resolve())

# ----------------------------------------------------------------------
# CLI/Dev Test Mode
# ----------------------------------------------------------------------
if __name__ == "__main__":
    builder = TemperatureBuilder(
        json_file="const_temp_param.json",
        fortran_input="ZDPlaskin/dbd0D.F90",
        fortran_output="ZDPlaskin/dbd0D.F90"
    )
    builder.update()
