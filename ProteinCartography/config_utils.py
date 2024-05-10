import enum
import pathlib


class ProteinCartographyInputError(Exception):
    pass


class Mode(enum.Enum):
    """
    The mode in which the pipeline is being run. There are two modes:
    1) 'search' mode: the pipeline will search for similar proteins to each of the input proteins
    using both blast and foldseek, then download structures from alphafold for the resulting hits,
    and finally cluster both the hits and the input proteins using foldseek
    2) 'cluster' mode: the pipeline will simply cluster all of the input proteins
    (for which the user must provide PDB files)
    """

    SEARCH = "search"
    CLUSTER = "cluster"


def _get_protids(config):
    """
    Determine the search-mode protein IDs and the 'key' protein IDs,
    given the mode in which the pipeline is being run and the user-provided config file
    """
    mode = Mode(config["mode"])
    input_dir = pathlib.Path(config["input_dir"])

    search_mode_input_protids = []
    key_protids = []

    if mode == Mode.SEARCH:
        # in search mode, the fasta files define the input protids
        # the list of input protein sequences (fasta files)
        input_fasta_filepaths = [
            filepath
            for filepath in input_dir.glob("*")
            if filepath.suffix[1:].lower() in ["fasta", "fa", "fna", "faa"]
        ]

        # there must be at least one such fasta file
        if not input_fasta_filepaths:
            raise ProteinCartographyInputError(
                "In 'search' mode, at least one FASTA file must be provided in the input directory."
            )

        search_mode_input_protids = [filepath.stem for filepath in input_fasta_filepaths]

        # in search mode, the 'key' protids are simply the input protids
        key_protids = search_mode_input_protids.copy()

    elif mode == Mode.CLUSTER:
        # in cluster mode, the only input files are PDB files (any fasta files are ignored)
        input_pdb_filepaths = [
            filepath for filepath in input_dir.glob("*") if filepath.suffix[1:].lower() == "pdb"
        ]

        # check that there is at least a reasonable number of PDB files provided
        # (enough that it makes sense to do the clustering)
        # TODO (KC): decide on a less arbitrary minimum number of PDBs
        min_num_pdb_files = 10
        if len(input_pdb_filepaths) < min_num_pdb_files:
            raise ProteinCartographyInputError(
                f"In 'cluster' mode, at least {min_num_pdb_files} PDB files must be provided "
                "in the input directory."
            )

        # in cluster mode, the 'key' protids are user-defined
        key_protids = config.get("key_protids", [])

        # check that the key protids are a subset of the protids for which PDB files were provided
        pdb_protids = [filepath.stem for filepath in input_pdb_filepaths]
        if not set(key_protids).issubset(pdb_protids):
            raise ProteinCartographyInputError(
                "The list of key proteins must be a subset of the list of input proteins."
            )

    return search_mode_input_protids, key_protids


def _get_features_file(config):
    """
    the features file is specific to 'cluster' mode and is required in that mode
    (because it is the source of the "Protein names" column required in `plot_semantic_analysis`)
    """
    mode = Mode(config["mode"])
    if mode != Mode.CLUSTER:
        return None

    features_file = config.get("features_file")
    if features_file is None:
        raise ProteinCartographyInputError(
            "A features file is required when running in cluster mode."
        )
    features_file = pathlib.Path(config["input_dir"]) / features_file
    if not features_file.is_file():
        raise ProteinCartographyInputError(
            f"No features file named '{features_file.name}' found in the input directory."
        )
    return features_file


def _get_features_override_file(config):
    """
    the override file is an optional user-provided file that can be used in either mode
    (see config.yml for more details)
    """
    features_override_file = pathlib.Path(config["input_dir"]) / config.get(
        "features_override_file", ""
    )

    # if the override file doesn't exist, ignore it
    # note: `features_override_file` cannot be set to `None`
    # because it is passed to a CLI option in the `aggregate_features` rule,
    # and snakemake serializes `None` to 'None';
    # instead, we must use an empty string, so that the CLI option is passed no value
    if not features_override_file.is_file():
        features_override_file = ""
    return features_override_file
