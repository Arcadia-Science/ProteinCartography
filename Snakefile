import os
from pathlib import Path


# Default pipeline configuration parameters are in this file
# If you create a new yml file and use the --configfile flag,
# options in that new file overwrite the defaults.
configfile: "./config.yml"


# Set the input directory
INPUT_DIR = Path(config["input_dir"])

# put most things into the output directory
OUTPUT_DIR = Path(config["output_dir"])

# Set the prefix of the output file of the analysis
ANALYSIS_NAME = config["analysis_name"]

# an optional user-provided override file
OVERRIDE_FILE = INPUT_DIR / config.get("override_file", "")

# if the override file doesn't exist, ignore it
# note: `OVERRIDE_FILE` cannot be set to `None` because it is passed to a CLI option
# in the `aggregate_features` rule, and snakemake serializes `None` to 'None';
# instead, we must use an empty string, so that the CLI option is passed no value
if not OVERRIDE_FILE.is_file():
    OVERRIDE_FILE = ""

TAXON_FOCUS = config.get("taxon_focus", "euk")

UNIPROT_ADDITIONAL_FIELDS = config["uniprot_additional_fields"]

MAX_BLAST_HITS = int(config["max_blast_hits"])
MAX_FOLDSEEK_HITS = int(config["max_foldseek_hits"])

BLAST_EVALUE = float(config["blast_evalue"])
BLAST_WORD_SIZE = int(config["blast_word_size"])
BLAST_WORD_SIZE_BACKOFF = int(config["blast_word_size_backoff"])
BLAST_NUM_ATTEMPTS = int(config["blast_num_attempts"])

MAX_STRUCTURES = int(config["max_structures"])

FS_DATABASES = config["foldseek_databases"]
FS_SERVER_ARG = f"--server {config['foldseek_server']}" if "foldseek_server" in config else ""
MODES = config["plotting_modes"]

MIN_LENGTH = int(config["min_length"])
MAX_LENGTH = int(config["max_length"])

BENCHMARKS_DIR = OUTPUT_DIR / "benchmarks"

# results from running blastp with the input proteins
BLAST_RESULTS_DIR = OUTPUT_DIR / "blast_results"

# results from calling the foldseek web API with the input proteins
FOLDSEEK_RESULTS_DIR = OUTPUT_DIR / "foldseek_results"

# metadata related to the hits (proteins) from blast and foldseek
PROTEIN_FEATURES_DIR = OUTPUT_DIR / "protein_features"

# PDBs downloaded from AlphaFold
PROTEIN_STRUCTURES_DIR = OUTPUT_DIR / "protein_structures"

# results from running foldseek to cluster the PDBs
FOLDSEEK_CLUSTERING_DIR = OUTPUT_DIR / "foldseek_clustering_results"

# final output results (plots and aggregated TSV files)
FINAL_RESULTS_DIR = OUTPUT_DIR / "final_results"

# create a list of protein IDs from the FASTA files in the input directory,
# allowing for multiple FASTA file formats
FASTA_FORMATS = [".fa", ".fasta", ".fna", ".faa"]
PROTIDS = []
for filename in os.listdir(INPUT_DIR):
    if any(filename.lower().endswith(suffix) for suffix in FASTA_FORMATS):
        protid = os.path.splitext(filename)[0]
        PROTIDS.append(protid)

BLAST_OUTPUT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "sacc",
    "saccver",
    "sgi",
    "staxids",
    "scomnames",
]
BLAST_OUTFMT = '"' + " ".join(["6"] + BLAST_OUTPUT_FIELDS) + '"'


wildcard_constraints:
    # constrain the `protid` wildcard to an ID in `PROTIDS`
    protid="|".join(PROTIDS),
    # constrain the `mode` wildcard to one of the plotting modes specified in the config file
    mode="|".join(MODES),


rule make_pdb:
    """
    Use the ESMFold API query to generate a pdb from a fasta file.
    This rule is only touched if a .pdb file doesn't already exist.
    """
    input:
        fasta_file=INPUT_DIR / "{protid}.fasta",
    output:
        pdb_file=INPUT_DIR / "{protid}.pdb",
    benchmark:
        BENCHMARKS_DIR / "{protid}.make_pdb.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/esmfold_apiquery.py -i {input.fasta_file}
        """


rule copy_pdb:
    """
    Copies existing or generated PDBs to the protein structures folder.
    """
    input:
        INPUT_DIR / "{protid}.pdb",
    output:
        PROTEIN_STRUCTURES_DIR / "{protid}.pdb",
    shell:
        """
        cp {input} {output}
        """


rule run_blast:
    """
    Using files located in the input directory, run `blastp` using the remote BLAST API.

    Large proteins will cause remote BLAST to fail; you can still perform a manual BLAST search to get around this.
    """
    input:
        fasta_file=INPUT_DIR / "{protid}.fasta",
    output:
        blast_results=BLAST_RESULTS_DIR / "{protid}.blast_results.tsv",
    benchmark:
        BENCHMARKS_DIR / "{protid}.run_blast.txt"
    conda:
        "envs/blast.yml"
    shell:
        """
        python ProteinCartography/run_blast.py \
          --query {input.fasta_file} \
          --out {output.blast_results} \
          --max_target_seqs {MAX_BLAST_HITS} \
          --outfmt {BLAST_OUTFMT} \
          --word_size {BLAST_WORD_SIZE} \
          --word_size_backoff {BLAST_WORD_SIZE_BACKOFF} \
          --num_attempts {BLAST_NUM_ATTEMPTS} \
          --evalue {BLAST_EVALUE}
        """


rule extract_blast_hits:
    """
    Extracts the top hits from the BLAST results file.
    """
    input:
        blast_results=BLAST_RESULTS_DIR / "{protid}.blast_results.tsv",
    output:
        blast_hits=BLAST_RESULTS_DIR / "{protid}.blast_hits.refseq.txt",
    benchmark:
        BENCHMARKS_DIR / "{protid}.extract_blast_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_blast_hits.py -i {input.blast_results} -o {output.blast_hits} -B {BLAST_OUTFMT}
        """


rule map_refseq_ids:
    """
    Map a list of RefSeq IDs to UniProt IDs using the Uniprot ID mapping API or bioservices.
    """
    input:
        blast_hits=rules.extract_blast_hits.output.blast_hits,
    output:
        blast_hits_uniprot_ids=BLAST_RESULTS_DIR / "{protid}.blast_hits.uniprot.txt",
    benchmark:
        BENCHMARKS_DIR / "{protid}.map_refseq_ids.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/map_refseq_ids.py -i {input.blast_hits} -o {output.blast_hits_uniprot_ids}
        """


rule run_foldseek:
    """
    Runs Foldseek using a query to the web API using a custom Python script.
    The script accepts an input file ending in '.pdb' and returns an output file ending in '.tar.gz'.
    The script also accepts a `--mode` flag of either '3diaa' (default) or 'tmalign' and choice of databases.
    After running, untars the files and extracts hits.

    Note: the foldseek web API returns a limited number of hits; up to 1000 per database
    """
    input:
        pdb_file=INPUT_DIR / "{protid}.pdb",
    output:
        foldseek_output=FOLDSEEK_RESULTS_DIR / "{protid}.fsresults.tar.gz",
        m8_files_dir=directory(FOLDSEEK_RESULTS_DIR / "{protid}"),
        m8_files=expand(FOLDSEEK_RESULTS_DIR / "{{protid}}" / "alis_{db}.m8", db=FS_DATABASES),
    params:
        fs_databases=expand("{fs_databases}", fs_databases=FS_DATABASES),
    conda:
        "envs/web_apis.yml"
    benchmark:
        BENCHMARKS_DIR / "{protid}.run_foldseek.txt"
    shell:
        """
        python ProteinCartography/foldseek_apiquery.py -i {input.pdb_file} -o {output.foldseek_output} {FS_SERVER_ARG} -d {params.fs_databases}
        tar -xvf {output.foldseek_output} -C {output.m8_files_dir}
        """


rule extract_foldseek_hits:
    input:
        m8_files=rules.run_foldseek.output.m8_files,
    output:
        foldseek_hits=FOLDSEEK_RESULTS_DIR / "{protid}.foldseek_hits.txt",
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_foldseek_hits.py -i {input.m8_files} -o {output.foldseek_hits} -m {MAX_FOLDSEEK_HITS}
        """


rule aggregate_foldseek_fraction_seq_identity:
    """
    Pulls the foldseek fraction sequence identity (fident) from the Foldseek results files for each input protein.

    This will probably be replaced in the future by an all-v-all sequence identity comparison using FAMSA, WITCH, or other approach.
    """
    input:
        m8_files=rules.run_foldseek.output.m8_files,
    output:
        fident_features=PROTEIN_FEATURES_DIR / "{protid}_fident_features.tsv",
    params:
        protid="{protid}",
    benchmark:
        BENCHMARKS_DIR / "{protid}.aggregate_foldseek_fident.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/aggregate_foldseek_fraction_seq_identity.py -i {input.m8_files} -o {output.fident_features} -p {params.protid}
        """


rule aggregate_hits:
    """
    Take all Uniprot ID lists and make them one big ID list, removing duplicates.
    """
    input:
        expand(rules.map_refseq_ids.output.blast_hits_uniprot_ids, protid=PROTIDS),
        expand(rules.extract_foldseek_hits.output.foldseek_hits, protid=PROTIDS),
    output:
        aggregated_hits=PROTEIN_FEATURES_DIR / "aggregated_hits.txt",
    benchmark:
        BENCHMARKS_DIR / "aggregate_hits.txt"
    shell:
        """
        python ProteinCartography/aggregate_hits.py -i {input} -o {output}
        """


rule fetch_uniprot_metadata:
    """
    Query Uniprot for the aggregated hits and download all metadata as a big ol' TSV.
    """
    input:
        rules.aggregate_hits.output.aggregated_hits,
    output:
        uniprot_features=PROTEIN_FEATURES_DIR / "uniprot_features.tsv",
    benchmark:
        BENCHMARKS_DIR / "fetch_uniprot_metadata.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py -i {input} -o {output.uniprot_features} -a {UNIPROT_ADDITIONAL_FIELDS}
        """


rule filter_aggregated_hits:
    """
    Use the metadata features from Uniprot to filter hits based on sequence status, fragment, and size.
    """
    input:
        rules.fetch_uniprot_metadata.output.uniprot_features,
    output:
        filtered_aggregated_hits=PROTEIN_FEATURES_DIR / "filtered_aggregated_hits.txt",
    benchmark:
        BENCHMARKS_DIR / "filter_aggregated_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/filter_aggregated_hits.py -i {input} -o {output.filtered_aggregated_hits} -m {MIN_LENGTH} -M {MAX_LENGTH}
        """


checkpoint download_pdbs:
    """
    Download all PDB files from AlphaFold
    """
    input:
        rules.filter_aggregated_hits.output.filtered_aggregated_hits,
    output:
        protein_structures_dir=directory(PROTEIN_STRUCTURES_DIR),
    benchmark:
        BENCHMARKS_DIR / "download_pdbs.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py -i {input} -o {output.protein_structures_dir} -M {MAX_STRUCTURES}
        """


def checkpoint_download_pdbs(wildcards):
    """
    Returns the paths to all PDB files
    (both those downloaded in `download_pdbs` and those corresponding to the input proteins)
    """
    # the directory containing the PDB files
    pdb_dirpath = checkpoints.download_pdbs.get(**wildcards).output.protein_structures_dir

    # the paths to the PDB files downloaded in `download_pdbs`
    pdb_filepaths = sorted(Path(pdb_dirpath).glob("*.pdb"))

    # append the paths to the PDB files corresponding to the input proteins
    # note: this triggers the `copy_pdb` rule to copy the input PDB files from `INPUT_DIR`
    pdb_filepaths += expand(PROTEIN_STRUCTURES_DIR / "{protid}.pdb", protid=PROTIDS)
    return pdb_filepaths


rule assess_pdbs:
    """
    Calculates the quality of all PDBs
    (those downloaded from AlphaFold and those corresponding to the input proteins).
    """
    input:
        checkpoint_download_pdbs,
    output:
        pdb_features=PROTEIN_FEATURES_DIR / "pdb_features.tsv",
    benchmark:
        BENCHMARKS_DIR / "assess_pdbs.txt"
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py -i {PROTEIN_STRUCTURES_DIR} -o {output.pdb_features}
        """


rule foldseek_clustering:
    """
    Runs foldseek all-v-all TM-score comparison and foldseek clustering.
    """
    input:
        checkpoint_download_pdbs,
    output:
        all_by_all_tmscores=FOLDSEEK_CLUSTERING_DIR / "all_by_all_tmscore_pivoted.tsv",
        struclusters_features=FOLDSEEK_CLUSTERING_DIR / "struclusters_features.tsv",
    conda:
        "envs/foldseek.yml"
    resources:
        mem_mb=32 * 1000,
    threads: 16
    benchmark:
        BENCHMARKS_DIR / "foldseek_clustering.txt"
    shell:
        """
        python ProteinCartography/foldseek_clustering.py -q {PROTEIN_STRUCTURES_DIR} -r {FOLDSEEK_CLUSTERING_DIR}
        """


rule dim_reduction:
    """
    Perform dimensionality reduction, saving as an embedding matrix and a TSV
    Write a set of functions to return Dataframes for interactive compute
    Write helper functions to save the dataframes only called by main()
    """
    input:
        rules.foldseek_clustering.output.all_by_all_tmscores,
    output:
        all_by_all_tmscores=FOLDSEEK_CLUSTERING_DIR / "all_by_all_tmscore_pivoted_{mode}.tsv",
    params:
        mode="{mode}",
    conda:
        "envs/analysis.yml"
    benchmark:
        BENCHMARKS_DIR / "{mode}.dim_reduction.txt"
    shell:
        """
        python ProteinCartography/dim_reduction.py -i {input} -m {params.mode}
        """


rule leiden_clustering:
    """
    Performs Leiden clustering on the data using scanpy's implementation.
    """
    input:
        rules.foldseek_clustering.output.all_by_all_tmscores,
    output:
        leiden_features=FOLDSEEK_CLUSTERING_DIR / "leiden_features.tsv",
    conda:
        "envs/analysis.yml"
    benchmark:
        BENCHMARKS_DIR / "leiden_clustering.txt"
    shell:
        """
        python ProteinCartography/leiden_clustering.py -i {input} -o {output.leiden_features}
        """


rule extract_input_protein_distances:
    """
    Extracts the distances from input proteins to other proteins in the dataset.
    """
    input:
        rules.foldseek_clustering.output.all_by_all_tmscores,
    output:
        distance_features=PROTEIN_FEATURES_DIR / "{protid}_distance_features.tsv",
    params:
        protid="{protid}",
    benchmark:
        BENCHMARKS_DIR / "{protid}.input_distances.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_input_protein_distances.py -i {input} -o {output.distance_features} -p {params.protid}
        """


rule calculate_concordance:
    """
    Currently, this subtracts the fraction sequence identity from the TM-score
    to get a measure of whether something is more similar in sequence or structure.

    We're working on developing some kind of test statistic that evaluates the significance
    of this difference from some expectation.
    """
    input:
        distance_features=rules.extract_input_protein_distances.output.distance_features,
        fident_features=rules.aggregate_foldseek_fraction_seq_identity.output.fident_features,
    params:
        protid="{protid}",
    output:
        concordance_features=PROTEIN_FEATURES_DIR / "{protid}_concordance_features.tsv",
    benchmark:
        BENCHMARKS_DIR / "{protid}.calculate_concordance.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/calculate_concordance.py -t {input.distance_features} -f {input.fident_features} -o {output.concordance_features} -p {params.protid}
        """


rule get_source_of_hits:
    """
    Checks the blast hits and foldseek hits to determine the source of each protein.
    """
    input:
        tm_scores=rules.foldseek_clustering.output.all_by_all_tmscores,
        hit_files=(
            expand(rules.extract_foldseek_hits.output.foldseek_hits, protid=PROTIDS)
            + expand(rules.map_refseq_ids.output.blast_hits_uniprot_ids, protid=PROTIDS)
        ),
    output:
        source_features=PROTEIN_FEATURES_DIR / "source_features.tsv",
    benchmark:
        BENCHMARKS_DIR / "get_source_of_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/get_source_of_hits.py -i {input.tm_scores} -f {input.hit_files} -o {output.source_features} -k {PROTIDS}
        """


rule aggregate_features:
    """
    Aggregate all TSV features provided by user in some specific directory, making one big TSV
    """
    input:
        rules.foldseek_clustering.output.struclusters_features,
        rules.fetch_uniprot_metadata.output.uniprot_features,
        rules.leiden_clustering.output.leiden_features,
        expand(rules.extract_input_protein_distances.output.distance_features, protid=PROTIDS),
        expand(
            rules.aggregate_foldseek_fraction_seq_identity.output.fident_features, protid=PROTIDS
        ),
        expand(rules.calculate_concordance.output.concordance_features, protid=PROTIDS),
        rules.get_source_of_hits.output.source_features,
        rules.assess_pdbs.output.pdb_features,
    output:
        aggregated_features=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_aggregated_features.tsv",
    benchmark:
        BENCHMARKS_DIR / "aggregate_features.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/aggregate_features.py -i {input} -o {output.aggregated_features} -v {OVERRIDE_FILE}
        """


rule plot_interactive:
    """
    Generate interactive scatter plot HTML programmatically based on user-input parameters
    Takes the TSV from rule aggregate_features and select default columns
    User should be able to call this module and pass their own functions to parse particular TSV columns
    Should have means to set a palette for each individual plot type, maybe as JSON?
    """
    input:
        tm_scores=rules.dim_reduction.output.all_by_all_tmscores,
        features=rules.aggregate_features.output.aggregated_features,
    output:
        html=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_aggregated_features_{{mode}}.html",
    params:
        mode="{mode}",
        protids=expand("{protid}", protid=PROTIDS),
    conda:
        "envs/plotting.yml"
    benchmark:
        BENCHMARKS_DIR / "{mode}.plot_interactive.txt"
    shell:
        """
        python ProteinCartography/plot_interactive.py -d {input.tm_scores} -f {input.features} -o {output.html} -t {params.mode} -k {params.protids} -x {TAXON_FOCUS}
        """


rule plot_similarity_leiden:
    """
    Plots a similarity score matrix for Leiden clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.
    """
    input:
        tm_scores=rules.foldseek_clustering.output.all_by_all_tmscores,
        features=rules.leiden_clustering.output.leiden_features,
    output:
        tsv=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_leiden_similarity.tsv",
        html=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_leiden_similarity.html",
    params:
        column="LeidenCluster",
    conda:
        "envs/plotting.yml"
    benchmark:
        BENCHMARKS_DIR / "plot_similarity_leiden.txt"
    shell:
        """
        python ProteinCartography/plot_cluster_similarity.py -m {input.tm_scores} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        """


rule plot_similarity_strucluster:
    """
    Plots a similarity score matrix for Foldseek's structural clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.

    TODO (KC): this rule is almost identical to `plot_similarity_leiden`
    """
    input:
        tm_scores=rules.foldseek_clustering.output.all_by_all_tmscores,
        features=rules.foldseek_clustering.output.struclusters_features,
    output:
        tsv=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_strucluster_similarity.tsv",
        html=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_strucluster_similarity.html",
    params:
        column="StruCluster",
    conda:
        "envs/plotting.yml"
    benchmark:
        BENCHMARKS_DIR / "plot_similarity_strucluster.txt"
    shell:
        """
        python ProteinCartography/plot_cluster_similarity.py -m {input.tm_scores} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        """


rule plot_semantic_analysis:
    """
    Plots a semantic analysis chart for groups within the data.
    """
    input:
        features=rules.aggregate_features.output.aggregated_features,
    output:
        pdf=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_semantic_analysis.pdf",
        html=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_semantic_analysis.html",
    params:
        agg_column="LeidenCluster",
        annot_column="'Protein names'",
    conda:
        "envs/plotting.yml"
    benchmark:
        BENCHMARKS_DIR / "plot_semantic_analysis.txt"
    shell:
        """
        python ProteinCartography/semantic_analysis.py -f {input.features} -c {params.agg_column} -n {params.annot_column} -o {output.pdf} -i {output.html} -a {ANALYSIS_NAME}
        """


rule plot_cluster_distributions:
    """
    Plots distributions of key values per cluster for each input protein.
    """
    input:
        features=rules.aggregate_features.output.aggregated_features,
    output:
        svg=FINAL_RESULTS_DIR / f"{ANALYSIS_NAME}_{{protid}}_distribution_analysis.svg",
    params:
        protid="{protid}",
    conda:
        "envs/plotting.yml"
    benchmark:
        BENCHMARKS_DIR / "plot_cluster_distributions_{protid}.txt"
    shell:
        """
        python ProteinCartography/plot_cluster_distributions.py -i {input.features} -o {output.svg} -k {params.protid}
        """


rule all:
    """
    This is a pseudo-rule that defines the final outputs of the pipeline

    Note: this rule appears at the end, rather than the beginning, of the snakefile
    in order to allow the definition of its inputs in terms of the outputs of other rules
    (whose definitions must appear before this rule)
    """
    # we use `default_target` to tell snakemake that this is the first rule to run
    # (it otherwise defaults to running the first rule in the snakefile)
    default_target: True
    input:
        rules.plot_similarity_leiden.output.html,
        rules.plot_similarity_strucluster.output.html,
        rules.plot_semantic_analysis.output.html,
        rules.plot_semantic_analysis.output.pdf,
        expand(rules.plot_interactive.output.html, mode=MODES),
        expand(rules.plot_cluster_distributions.output.svg, protid=PROTIDS),
