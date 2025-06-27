#!/usr/bin/env python3
"""
Filter Kinet File Based on Selected Reactions

This script constructs a reduced version of `kinet.inp` based on a list of
approved reactions (`reacciones_usadas.txt`), preserving:
- Direct reactions
- Grouped reactions (with @ clusters)
- BOLSIG+ reactions with clustering
"""

import re
import argparse
from pathlib import Path
import logging

LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)

# ----------------------------------------------------------------------
# Load Reactions to Keep
# ----------------------------------------------------------------------
def cargar_reacciones_objetivo(path_txt):
    """
    Load the reduced reaction list from a text file.

    Returns:
        reacciones (set): Standard reactions (as strings).
        bolsig_reacciones (set): BOLSIG+ reactions (special handling).
    """
    with open(path_txt, 'r', encoding='utf-8') as file:
        raw_lines = file.readlines()

    reacciones = set()
    bolsig_reacciones = set()

    for line in raw_lines:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        if "bolsig:" in parts[1].lower():
            bolsig_reacciones.add(parts[1].lower().replace(" ", ""))
        else:
            reaccion = " ".join(parts[1:]).replace(" ", "").lower()
            reacciones.add(reaccion)

    return reacciones, bolsig_reacciones

# ----------------------------------------------------------------------
# Group Handling: Extract and Replace Clusters
# ----------------------------------------------------------------------
def extraer_clusters_posicionales(lines, start_idx):
    """
    Extracts cluster (grouped species) definitions starting at a given line.

    Returns:
        clusters_raw (dict): e.g. {'@X': ['NO', 'NO2']}
        clusters_lower (dict): lower-case version
        next_index (int): next index to resume parsing
    """
    clusters_raw = {}
    clusters_lower = {}
    i = start_idx

    while i < len(lines):
        line = lines[i]
        if not line.strip().startswith("@"):
            break
        nombre, valores = line.strip().split("=", 1)
        nombre = nombre.strip()
        valores_raw = valores.strip().split()
        clusters_raw[nombre] = valores_raw
        clusters_lower[nombre] = [v.lower() for v in valores_raw]
        i += 1

    return clusters_raw, clusters_lower, i

def encontrar_at_vars(texto):
    """
    Find all @variables used in a string.
    """
    return re.findall(r"@[\w]+", texto)

# ----------------------------------------------------------------------
# Core Filtering Logic
# ----------------------------------------------------------------------
def filtrar_kinet(kinet_path, reacciones_final_path, salida_path):
    """
    Main function to filter a kinet.inp file based on selected reactions.

    Args:
        kinet_path (str): Path to input kinet.inp
        reacciones_final_path (str): Path to file with reactions to preserve
        salida_path (str): Output path for filtered kinet file
    """
    with open(kinet_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    reacciones_objetivo, reacciones_bolsig = cargar_reacciones_objetivo(reacciones_final_path)

    en_reactions = False
    nueva_kinet = []
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Handle outside/inside REACTIONS block
        if stripped.upper() == "REACTIONS":
            en_reactions = True
            nueva_kinet.append(line)
            i += 1
            continue

        if en_reactions and stripped.upper() == "END":
            en_reactions = False
            nueva_kinet.append(line)
            i += 1
            continue

        if not en_reactions:
            nueva_kinet.append(line)
            i += 1
            continue

        # Comments or directives — keep as-is
        if stripped.startswith("#") or stripped == "" or stripped.startswith("$"):
            nueva_kinet.append(line)
            i += 1
            continue

        # ------------------------------------------------------------------
        # Case 1: BOLSIG+ reactions (can be grouped with @)
        # ------------------------------------------------------------------
        elif "BOLSIG" in line.upper():
            match = re.match(r"(.+)=>(.+)!(.+)", line.strip())
            if not match:
                i += 1
                continue

            lhs_raw = match.group(1).strip()
            rhs_raw = match.group(2).strip()
            tag_raw = match.group(3).strip()

            bolsig_match = re.search(r"BOLSIG\s+(\S+)\s*->\s*(\S+)", tag_raw, re.IGNORECASE)
            if not bolsig_match:
                i += 1
                continue

            reactivo = bolsig_match.group(1).strip()
            producto = bolsig_match.group(2).strip()

            # Case: Simple reaction (no @)
            if not reactivo.startswith("@") and not producto.startswith("@"):
                clave = f"bolsig:{reactivo}->{producto}".lower().replace(" ", "")
                if clave in reacciones_bolsig:
                    nueva_kinet.append(line)
                i += 1
                continue

            # Case: Grouped reactions
            clusters_raw, clusters_lower, next_idx = extraer_clusters_posicionales(lines, i + 1)

            if producto.startswith("@"):
                grupo = clusters_lower.get(producto, [])
                aux = {k: v for k, v in clusters_raw.items() if k != producto}
                usadas = [idx for idx in range(len(grupo))
                          if f"bolsig:{reactivo}->{grupo[idx]}".lower().replace(" ", "") in reacciones_bolsig]

                if usadas:
                    nueva_kinet.append(line)
                    nueva_kinet.append(f"    {producto} = {'  '.join([clusters_raw[producto][i] for i in usadas])}\n")
                    for var, lista in aux.items():
                        valores = [lista[i] for i in usadas if i < len(lista)]
                        if valores:
                            nueva_kinet.append(f"    {var} = {'  '.join(valores)}\n")

                i = next_idx
                continue

            elif reactivo.startswith("@"):
                grupo = clusters_lower.get(reactivo, [])
                aux = {k: v for k, v in clusters_raw.items() if k != reactivo}
                usadas = [idx for idx in range(len(grupo))
                          if f"bolsig:{grupo[idx]}->{producto}".lower().replace(" ", "") in reacciones_bolsig]

                if usadas:
                    nueva_kinet.append(line)
                    nueva_kinet.append(f"    {reactivo} = {'  '.join([clusters_raw[reactivo][i] for i in usadas])}\n")
                    for var, lista in aux.items():
                        valores = [lista[i] for i in usadas if i < len(lista)]
                        if valores:
                            nueva_kinet.append(f"    {var} = {'  '.join(valores)}\n")

                i = next_idx
                continue

        # ------------------------------------------------------------------
        # Case 2: Grouped reactions with @ variables
        # ------------------------------------------------------------------
        elif "@" in line:
            match = re.match(r"(.+)=>(.+)!(.+)", line.strip())
            if not match:
                i += 1
                continue

            lhs_raw, rhs_raw, tag_raw = match.groups()
            clusters_raw, clusters_lower, next_idx = extraer_clusters_posicionales(lines, i + 1)

            variables_lhs = encontrar_at_vars(lhs_raw)
            variables_rhs = encontrar_at_vars(rhs_raw)
            variables_tag = encontrar_at_vars(tag_raw)

            all_vars = list(set(variables_lhs + variables_rhs))
            especies_base = clusters_lower.get(all_vars[0], []) if all_vars else []
            valid_indices = []

            for idx in range(len(especies_base)):
                lhs_i, rhs_i, tag_i = lhs_raw, rhs_raw, tag_raw
                valid = True

                for var in all_vars:
                    if idx >= len(clusters_lower[var]):
                        valid = False
                        break
                    lhs_i = lhs_i.replace(var, clusters_lower[var][idx])
                    rhs_i = rhs_i.replace(var, clusters_lower[var][idx])

                for var in variables_tag:
                    if idx >= len(clusters_raw[var]):
                        valid = False
                        break
                    tag_i = tag_i.replace(var, clusters_raw[var][idx])

                if not valid:
                    continue

                reac_norm = f"{lhs_i.strip()} => {rhs_i.strip()}".replace(" ", "").lower()
                if reac_norm in reacciones_objetivo:
                    valid_indices.append(idx)

            if valid_indices:
                nueva_kinet.append(line)
                for var in all_vars + variables_tag:
                    valores = [clusters_raw[var][j] for j in valid_indices if j < len(clusters_raw[var])]
                    if valores:
                        nueva_kinet.append(f"    {var} = {'  '.join(valores)}\n")

            i = next_idx
            continue

        # ------------------------------------------------------------------
        # Case 3: Simple reactions
        # ------------------------------------------------------------------
        elif "=>" in line and "!" in line:
            reaccion_izq = line.split("!")[0].replace(" ", "").lower()
            if reaccion_izq in reacciones_objetivo:
                nueva_kinet.append(line)
            i += 1
            continue

        # Fallback
        i += 1

    # Save output file
    with open(salida_path, 'w', encoding='utf-8') as f_out:
        f_out.writelines(nueva_kinet)

    logging.info(f"New kinet.inp file saved in: {salida_path}")


# ----------------------------------------------------------------------
# Entry Point for External Calls
# ----------------------------------------------------------------------
def run_new_kinet(
    
    kinet_path="Data/Reduced_Sets/kinet_original.inp",
    reacciones_path="reacciones_usadas.txt",
    output_path="kinet.inp"
):
    """
    Wrapper function to run kinet filtering with default paths.
    """
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    HELPER_FILES = PROJECT_ROOT / "src/helper_files"

    filtrar_kinet(PROJECT_ROOT / kinet_path, HELPER_FILES / reacciones_path, PROJECT_ROOT / "output" / output_path)
