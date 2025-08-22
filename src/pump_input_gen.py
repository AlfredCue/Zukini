#!/usr/bin/env python3
"""
PumpKin Input File Generator

This script creates an `input.txt` file compatible with PumpKin for a given species
by retrieving its corresponding ID from a QtPlasKin species list.
"""

import os
import sys

# ----------------------------------------------------------------------
# Helper: Get Species ID from QtPlasKin List
# ----------------------------------------------------------------------
def get_species_id(species_name: str, qt_species_list_path: str) -> str:
    """
    Finds the ID tag of the given species in the QtPlasKin species list.

    Args:
        species_name (str): The species to search for.
        qt_species_list_path (str): Path to the QtPlasKin species list file.

    Returns:
        str: The tag/ID associated with the species.

    Raises:
        FileNotFoundError: If the file is missing.
        ValueError: If the species name isn't found.
    """
    if not os.path.exists(qt_species_list_path):
        raise FileNotFoundError(f"File not found: {qt_species_list_path}")

    with open(qt_species_list_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                tag, name = parts[0], parts[1]
                if name.lower() == species_name.lower():
                    return tag

    raise ValueError(f"Species '{species_name}' not found in Qt species list.")

# ----------------------------------------------------------------------
# Helper: Write Input File
# ----------------------------------------------------------------------
def create_pumpkin_input_file(tag: str, species_name: str, output_dir: str = "pumpkin-1.4/src/Examples/Sim") -> None:
    """
    Writes a PumpKin-compatible input file for the specified species.

    Args:
        tag (str): Species tag ID.
        species_name (str): Species name (used in output filename).
        output_dir (str): Output directory for the input file.
    """
    os.makedirs(output_dir, exist_ok=True)

    input_template = (
        "interest \t\t\t= {interest}\n"
        "t_init \t\t\t\t= 4.8677e-06\n"
        "t_end \t\t\t\t= 1.799876E-01\n"
        "max_bp \t\t\t\t= 0\n"
        "tau_lifetime \t\t= 0.9e-4\n"
        "max_path \t\t\t= 500\n"
        "f_min \t\t\t\t= 0.01\n"
    )

    filename = os.path.join(output_dir, "input.txt")
    with open(filename, 'w') as f:
        f.write(input_template.format(interest=tag))

    print(f"[INFO] Created input file: {filename}")

# ----------------------------------------------------------------------
# Main Callable Function
# ----------------------------------------------------------------------
def generate_pumpkin_input(
    species: str,
    qt_species_list_path: str = "Data/QtFiles/Qt_1/qt_species_list.txt",
    output_dir: str = "pumpkin-1.4/src/Examples/Sim"
) -> None:
    """
    Main routine to generate a PumpKin input for a given species.

    Args:
        species (str): Target species name (case-insensitive).
        qt_species_list_path (str): QtPlasKin species list.
        output_dir (str): Directory where the input.txt file will be placed.
    """
    try:
        tag = get_species_id(species, qt_species_list_path)
        create_pumpkin_input_file(tag, species, output_dir)
    except Exception as e:
        print(f"[ERROR] {e}")
        raise

# ----------------------------------------------------------------------
# CLI Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("[ERROR] Usage: python Pump_input_gen.py <species_name>")
        sys.exit(1)

    species = sys.argv[1]
    generate_pumpkin_input(species)
