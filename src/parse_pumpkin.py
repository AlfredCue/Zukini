#!/usr/bin/env python3
"""
PumpKin Output Parser

Parses pathway reactions from PumpKin output, tags them in the master reaction list,
and extracts all species involved (excluding placeholders). Results are stored in:
- reactions_tagged.json (for filtering)
- species_to_expand.json (for next generation expansion)
"""

import json
import os
import re

PLACEHOLDER_SPECIES = {"any_neutral", "any_species", "any_ion_positive", "any_ion_negative"}

# ----------------------------------------------------------------------
# STEP 1: Extract Pathway Reactions from PumpKin Output
# ----------------------------------------------------------------------
def parse_pumpkin_output(filename):
    """
    Parses reaction lines from PumpKin's output text file.

    Returns:
        List of reaction strings extracted from pathways.
    """
    reactions = set()
    inside_pathway_table = False

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not inside_pathway_table and "|" in line and "Pathway" in line:
                inside_pathway_table = True
                continue
            if inside_pathway_table and line.startswith("################"):
                break
            if inside_pathway_table and line.startswith("|") and "(" in line and ")" in line:
                try:
                    start = line.find("(")
                    end = line.rfind(")")
                    reaction = line[start + 1:end].strip()
                    reactions.add(reaction)
                except Exception as e:
                    print(f"[ERROR] Failed to parse line: {line} — {e}")

    print(f"[INFO] Extracted {len(reactions)} unique reactions from pathways.")
    return list(reactions)

# ----------------------------------------------------------------------
# STEP 2: Load All Reactions from Qt Files
# ----------------------------------------------------------------------
def load_reaction_list(filepath):
    """
    Load all reactions from qt_reactions_list.txt (ignores indices).

    Returns:
        List of reaction strings.
    """
    reactions = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                reactions.append(parts[1].strip())
    return reactions

# ----------------------------------------------------------------------
# STEP 3: Tag Reactions (pathway-preserving)
# ----------------------------------------------------------------------
def tag_reactions(all_reactions, pathway_reactions, previous_tags=None):
    """
    Tag each reaction as True (if present in pathway) or preserve True from previous run.

    Args:
        all_reactions: list of all available reactions.
        pathway_reactions: list of reactions from PumpKin pathways.
        previous_tags: optional existing tags (dict).

    Returns:
        dict: {reaction_str: bool}
    """
    pathway_set = set(pathway_reactions)
    if previous_tags is None:
        previous_tags = {}

    tagged = {}
    for reaction in all_reactions:
        was_true = previous_tags.get(reaction, False)
        is_true_now = reaction in pathway_set
        tagged[reaction] = was_true or is_true_now

    return tagged

# ----------------------------------------------------------------------
# STEP 4: Extract All Unique Species (excluding placeholders)
# ----------------------------------------------------------------------
def extract_species_from_reactions(reaction_list):
    """
    Extract unique species (ignoring placeholder species and BOLSIGs).

    Returns:
        Set of lowercase species names.
    """
    species_set = set()

    for reaction in reaction_list:
        reaction = reaction.lower()
        if reaction.startswith("bolsig:"):
            continue
        if "=>" not in reaction:
            continue

        reactants, products = reaction.split("=>", 1)
        for side in (reactants, products):
            species = smart_species_split(side)
            for specie in species:
                specie = specie.strip()
                if specie and specie not in PLACEHOLDER_SPECIES:
                    species_set.add(specie)

    print(f"[INFO] Extracted {len(species_set)} unique species from reactions.")
    return species_set

# ----------------------------------------------------------------------
# Utility: Handle Ion Species with + or -
# ----------------------------------------------------------------------
def smart_species_split(side_str):
    """
    Splits a reaction side (e.g. 'NO+O2') into species, handling ^+ or ^- ions safely.

    Args:
        side_str (str): reactants or products.

    Returns:
        List of species strings.
    """
    species = []
    buffer = ""
    i = 0
    while i < len(side_str):
        if side_str[i] == '+':
            if i > 0 and side_str[i - 1] == '^':
                buffer += '+'
            else:
                if buffer:
                    species.append(buffer)
                buffer = ""
        else:
            buffer += side_str[i]
        i += 1
    if buffer:
        species.append(buffer)
    return species

# ----------------------------------------------------------------------
# Utility: Save JSON
# ----------------------------------------------------------------------
def save_as_json(data, filename="reactions_tagged.json"):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"[INFO] Saved tagged reactions to {filename}")

# ----------------------------------------------------------------------
# Master Runner
# ----------------------------------------------------------------------
def run_parse_pumpkin(
    qt_reactions_path="Data/QtFiles/Qt_1/qt_reactions_list.txt",
    pumpkin_output_path="pumpkin-1.4/pumpkin_output.txt",
    json_output_path="reactions_tagged.json",
    species_output_path="species_to_expand.json"
):
    """
    Main function that parses PumpKin pathways, updates reaction tags,
    and extracts species for next expansion step.
    """
    # Step 1: Extract pathway reactions
    pathway_reactions = parse_pumpkin_output(pumpkin_output_path)

    # Step 2: Load all available reactions
    all_reactions = load_reaction_list(qt_reactions_path)

    # Step 3: Load previous tags (if available)
    previous_tags = {}
    if os.path.exists(json_output_path):
        with open(json_output_path, 'r') as f:
            previous_tags = json.load(f)

    # Step 4: Update and save tags
    tagged_data = tag_reactions(all_reactions, pathway_reactions, previous_tags)
    save_as_json(tagged_data, json_output_path)

    # Step 5: Extract and save species to expand
    species_set = extract_species_from_reactions(pathway_reactions)
    with open(species_output_path, 'w') as f:
        json.dump(list(species_set), f, indent=4)

    print(f"[INFO] Saved extracted species to {species_output_path}")

# ----------------------------------------------------------------------
# CLI Entry Point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    run_parse_pumpkin()
