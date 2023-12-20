#!/usr/bin/env python
import argparse
import os
import re
from io import StringIO
from pathlib import Path

import arcadia_pycolor as apc
import numpy as np
import pandas as pd

__all__ = [
    "fetch_atoms",
    "fetch_dbref",
    "fetch_experiment",
    "fetch_remark",
    "fetch_title",
    "extract_residue_confidence",
    "assign_residue_colors",
    "parse_chains",
    "assign_origin",
    "assess_pdbs",
]

RESIDUE_CONFIDENCE_COLORS = {
    "very_high": "#4A72B0",
    "confident": apc.All["arcadia:vitalblue"],
    "low": apc.All["arcadia:canary"],
    "very_low": apc.All["arcadia:amber"],
}

RESIDUE_BINS = ["very_low", "low", "confident", "very_high"]

# based on: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
ATOM_SPEC_DICT = {
    "ATOM": (0, 6),
    "SERIAL": (6, 11),
    "NAME": (12, 16),
    "ALTLOC": (16, 17),
    "RESIDUE": (17, 20),
    "CHAIN": (21, 22),
    "RESNUM": (22, 26),
    "RESINS": (26, 27),
    "X": (30, 38),
    "Y": (38, 45),
    "Z": (46, 54),
    "OCC": (54, 60),
    "TEMP": (60, 66),
    "SEG": (72, 76),
    "ELEM": (76, 78),
    "CHARGE": (78, 80),
}

ATOM_SPEC_VALS = list(ATOM_SPEC_DICT.values())
ATOM_SPEC_NAMES = list(ATOM_SPEC_DICT.keys())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the directory of PDB files to assess",
    )
    parser.add_argument("-o", "--output", required=True, help="Name of output TSV file.")
    args = parser.parse_args()
    return args


def is_valid_pdb(input_path: str) -> bool:
    """
    Checks if a file is a valid PDB file.
    Note: if the PDB file was downloaded from Alphafold, this usually means that the protein
    does not have a structure in AlphaFold.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        return "<Error>" not in f.read()


def fetch_atoms(input_path: str) -> pd.DataFrame:
    """
    Retrieves atoms as a dataframe from a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        atoms = [i for i in f.readlines() if "ATOM" in i[0:6]]

    if len(atoms) == 0:
        return pd.DataFrame()

    data = pd.read_fwf(StringIO("".join(atoms)), names=ATOM_SPEC_NAMES, colspecs=ATOM_SPEC_VALS)

    return data


def fetch_dbref(input_path: str) -> pd.DataFrame:
    """
    Retrieves database references from a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        dbref = [i for i in f.readlines() if "DBREF" in i[0:6]]

    if len(dbref) == 0:
        return pd.DataFrame()

    data = pd.read_fwf(StringIO("".join(dbref)), header=None)

    return data


def fetch_experiment(input_path: str) -> str:
    """
    Retrieves experimental data from a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        expdta = [
            " ".join([i for i in i.split() if i != "EXPDTA"])
            for i in f.readlines()
            if "EXPDTA" in i
        ]

    expdta_out = " ".join(expdta).split(";")

    return expdta_out


def fetch_title(input_path: str) -> str:
    """
    Retrieves title of a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        title = [
            " ".join([i for i in i.split() if i != "TITLE"]) for i in f.readlines() if "TITLE" in i
        ]

    title_out = " ".join(title)

    return title_out


def fetch_remark(input_path: str) -> str:
    """
    Retrieves remarks from a PDB file and returns the contents as a concatenated string.

    Args:
        input_path (str): path of PDB file.
    """
    with open(input_path) as f:
        remark = [
            " ".join([i for i in i.split() if i != "REMARK"])
            for i in f.readlines()
            if "REMARK" in i
        ]

    remark_out = " ".join(remark)

    return remark_out


def extract_residue_confidence(input_path: str):
    """
    Extracts residue confidence ("TEMP" or temperature) from a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    data = fetch_atoms(input_path)

    return list(data["TEMP"].astype(float))


def assign_residue_colors(lst: list):
    """
    Assigns a color to each residue in a PDB based on the AlphaFold color scheme
    (converted to analogous Arcadia colors).

    Args:
        lst (list): list of atoms temperatures.
    """
    bins = [0, 50, 70, 90, 100]
    result = np.digitize(lst, bins, right=True)
    colors = [RESIDUE_CONFIDENCE_COLORS[RESIDUE_BINS[i - 1]] for i in result]

    return colors


def parse_chains(input_path: str):
    """
    Extracts list of chains from a PDB file.

    Args:
        input_path (str): path of PDB file.
    """
    data = fetch_atoms(input_path)
    chains = list(data["CHAIN"].unique())

    return chains


def assign_origin(input_path: str):
    """
    Assigns an origin to a PDB file based on presence/ absence of references to AlphaFold,
    Protein Data Bank, or ESMFold.

    Args:
        input_path (str): path of PDB file.
    """
    AF_FLAG, AF_TITLE_FLAG, AF_REMARK_FLAG = 0, 0, 0
    PDB_FLAG, PDB_REF_FLAG, PDB_REMARK_FLAG = 0, 0, 0
    ESM_FLAG, ESM_TITLE_FLAG, ESM_REMARK_FLAG = 0, 0, 0

    with open(input_path) as f:
        contents = f.read()

        if re.search("ALPHAFOLD", contents, re.IGNORECASE):
            AF_FLAG = 1
        if re.search("PDB", contents, re.IGNORECASE):
            PDB_FLAG = 1
        if re.search("ESMFOLD", contents, re.IGNORECASE):
            ESM_FLAG = 1
    try:
        dbref = fetch_dbref(input_path)[5].values
    except KeyError:
        dbref = []

    if "PDB" in dbref:
        PDB_REF_FLAG = 1

    title = fetch_title(input_path).upper()
    if "ALPHAFOLD" in title:
        AF_TITLE_FLAG = 1
    elif "ESMFOLD" in title:
        ESM_TITLE_FLAG = 1

    remark = fetch_remark(input_path).upper()
    if "ALPHAFOLD" in remark:
        AF_REMARK_FLAG = 1
    elif "ESMFOLD" in remark:
        ESM_REMARK_FLAG = 1
    elif "RCSB" in remark:
        PDB_REMARK_FLAG = 1

    AF_SCORE = AF_FLAG + AF_TITLE_FLAG + AF_REMARK_FLAG
    PDB_SCORE = PDB_FLAG + PDB_REF_FLAG + PDB_REMARK_FLAG
    ESM_SCORE = ESM_FLAG + ESM_TITLE_FLAG + ESM_REMARK_FLAG

    OTHER_SCORE = 3 - AF_SCORE - PDB_SCORE - ESM_SCORE

    scores = {
        "AlphaFold": AF_SCORE,
        "ESMFold": ESM_SCORE,
        "PDB": PDB_SCORE,
        "Other": OTHER_SCORE,
    }

    maxscore = max(zip(scores.values(), scores.keys()))[1]

    return maxscore


def assess_pdbs(structure_filepaths: list, output_file=None):
    """
    Assesses PDB quality, experimental information, origin,
    and lists chains for a list of PDB paths.

    Args:
        structure_filepaths (list): list of paths to the PDB files to assess.
        output_file (str): path to output file, if saving results.
    """
    collector_df = pd.DataFrame()

    for structure_filepath in structure_filepaths:
        if not os.path.exists(structure_filepath):
            continue

        if not is_valid_pdb(structure_filepath):
            continue

        origin = assign_origin(structure_filepath)

        if origin != "PDB":
            max_confidence = np.max(extract_residue_confidence(structure_filepath))
            min_confidence = np.min(extract_residue_confidence(structure_filepath))
            if max_confidence <= 1 and min_confidence <= 1:
                confidence = np.average(extract_residue_confidence(structure_filepath)) * 100
            elif max_confidence >= 1 and min_confidence >= 1:
                confidence = np.average(extract_residue_confidence(structure_filepath))
        else:
            confidence = 100

        chains = parse_chains(structure_filepath)

        info_df = pd.DataFrame(
            {
                "protid": os.path.basename(structure_filepath).split(".pdb")[0],
                "pdb_origin": [origin],
                "pdb_confidence": [confidence],
                "pdb_chains": [chains],
            }
        )

        collector_df = pd.concat([collector_df, info_df], ignore_index=True)

    if output_file is not None:
        collector_df.to_csv(output_file, sep="\t", index=None)

    return collector_df


def main():
    args = parse_args()
    structure_filepaths = Path(args.input).glob("*.pdb")
    assess_pdbs(structure_filepaths, output_file=args.output)


if __name__ == "__main__":
    main()
