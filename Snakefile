import os
from pathlib import Path


# Default pipeline configuration parameters are in this file
# If you create a new yml file and use the --configfile flag,
# options in that new file overwrite the defaults.
configfile: "./config.yml"


# Set the input directory
#### In the future, also accept Uniprot accession numbers, which will be auto-queried and downloaded
input_dir = Path(config["input_dir"])

# Set the prefix of the output file of the analysis
analysis_name = config["analysis_name"]

# put most things into the output directory
output_dir = Path(config["output_dir"])

# Check for an override file, setting a variable if it exists
# note: `OVERRIDE_FILE` cannot be `None` because it is passed to a CLI option
# in the `aggregate_features` rule, and snakemake serializes `None` to 'None';
# using an empty string instead results in the CLI option being passed no value
if "override_file" in config:
    OVERRIDE_FILE = input_dir / config["override_file"]

    # If it isn't a real file, ignore it
    if not os.path.exists(OVERRIDE_FILE):
        OVERRIDE_FILE = ""
else:
    OVERRIDE_FILE = ""

if "taxon_focus" in config:
    TAXON_FOCUS = config["taxon_focus"]
else:
    TAXON_FOCUS = "euk"

UNIPROT_ADDITIONAL_FIELDS = config["uniprot_additional_fields"]

MAX_BLASTHITS = int(config["max_blasthits"])
MAX_FOLDSEEKHITS = int(config["max_foldseekhits"])

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

# these directories fall within the output directory
blastresults_dir = Path("blastresults/")
foldseekresults_dir = Path("foldseekresults/")
foldseekclustering_dir = Path("foldseekclustering/")

# directory in which to copy or download PDB files
pdb_download_dir = output_dir / foldseekclustering_dir

clusteringresults_dir = Path("clusteringresults/")
benchmarks_dir = Path("benchmarks/")

# gets the protein ID based on FASTA file name
# flexibly checks if fasta file is correct suffix
FASTA_FORMATS = [".fa", ".fasta", ".fna", ".faa"]
PROTID = []
for file in os.listdir(input_dir):
    if any(file.lower().endswith(suffix) for suffix in FASTA_FORMATS):
        file_id = os.path.splitext(file)[0]
        PROTID.append(file_id)

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


rule all:
    input:
        expand(
            output_dir
            / clusteringresults_dir
            / (analysis_name + "_aggregated_features_{modes}.html"),
            modes=MODES,
        ),
        output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.html"),
        output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.html"),
        output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.pdf"),
        output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.html"),
        expand(
            output_dir
            / clusteringresults_dir
            / (analysis_name + "_{protid}_distribution_analysis.svg"),
            protid=PROTID,
        ),


rule make_pdb:
    """
    Use the ESMFold API query to generate a pdb from a fasta file.
    This rule is only touched if a .pdb file doesn't already exist.
    """
    input:
        cds=input_dir / "{protid}.fasta",
    output:
        pdb=pdb_download_dir / "{protid}.pdb",
    benchmark:
        output_dir / benchmarks_dir / "{protid}.make_pdb.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/esmfold_apiquery.py -i {input.cds}
        """


rule copy_pdb:
    """
    Copies existing or generated PDBs to the Foldseek clustering folder.
    """
    input:
        input_dir / "{protid}.pdb",
    output:
        pdb_download_dir / "{protid}.pdb",
    shell:
        """
        cp {input} {output}
        """


# first try to copy any user-provided PDB files from the input directory;
# if they don't exist, generate them using make_pdb
ruleorder: copy_pdb > make_pdb


rule run_blast:
    """
    Using files located in the input directory, run `blastp` using the remote BLAST API.

    Large proteins will cause remote BLAST to fail; you can still perform a manual BLAST search to get around this.
    """
    input:
        fasta_file=input_dir / "{protid}.fasta",
    output:
        blast_results=output_dir / blastresults_dir / "{protid}.blastresults.tsv",
    params:
        max_target_seqs=MAX_BLASTHITS,
        outfmt=BLAST_OUTFMT,
        word_size=BLAST_WORD_SIZE,
        word_size_backoff=BLAST_WORD_SIZE_BACKOFF,
        evalue=BLAST_EVALUE,
        num_attempts=BLAST_NUM_ATTEMPTS,
    benchmark:
        output_dir / benchmarks_dir / "{protid}.run_blast.txt"
    conda:
        "envs/blast.yml"
    shell:
        """
        python ProteinCartography/run_blast.py \
          --query {input.fasta_file} \
          --out {output.blast_results} \
          --max_target_seqs {params.max_target_seqs} \
          --outfmt {params.outfmt} \
          --word_size {params.word_size} \
          --word_size_backoff {params.word_size_backoff} \
          --num_attempts {params.num_attempts} \
          --evalue {params.evalue}
        """


rule extract_blast_hits:
    """
    Extracts the top hits from the BLAST results file.
    """
    input:
        blastresults=output_dir / blastresults_dir / "{protid}.blastresults.tsv",
    output:
        refseqhits=output_dir / blastresults_dir / "{protid}.blasthits.refseq.txt",
    params:
        blast_outfmt=BLAST_OUTFMT,
    benchmark:
        output_dir / benchmarks_dir / "{protid}.extract_blast_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_blasthits.py -i {input.blastresults} -o {output.refseqhits} -B {params.blast_outfmt}
        """


rule map_refseqids:
    """
    Using List of RefSeq IDs, query the Uniprot ID mapping tool using bioservices UniProt.
    Returns a list of UniProt IDs.
    """
    input:
        refseqhits=output_dir / blastresults_dir / "{protid}.blasthits.refseq.txt",
    output:
        uniprothits=output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt",
    benchmark:
        output_dir / benchmarks_dir / "{protid}.map_refseqids.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/map_refseqids.py -i {input.refseqhits} -o {output.uniprothits}
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
        cds=input_dir / "{protid}.pdb",
    output:
        targz=output_dir / foldseekresults_dir / "{protid}.fsresults.tar.gz",
        unpacked=directory(output_dir / foldseekresults_dir / "{protid}"),
        m8files=expand(
            output_dir / foldseekresults_dir / "{{protid}}" / "alis_{db}.m8", db=FS_DATABASES
        ),
    params:
        fs_databases=expand("{fs_databases}", fs_databases=FS_DATABASES),
        fs_server_arg=FS_SERVER_ARG,
    conda:
        "envs/web_apis.yml"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.run_foldseek.txt"
    shell:
        """
        python ProteinCartography/foldseek_apiquery.py -i {input.cds} -o {output.targz} {params.fs_server_arg} -d {params.fs_databases}
        tar -xvf {output.targz} -C {output.unpacked}
        """


rule extract_foldseek_hits:
    input:
        m8files=rules.run_foldseek.output.m8files,
    output:
        foldseekhits=output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt",
    params:
        max_foldseekhits=MAX_FOLDSEEKHITS,
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_foldseekhits.py -i {input.m8files} -o {output.foldseekhits} -m {params.max_foldseekhits}
        """


rule aggregate_foldseek_fraction_seq_identity:
    """
    Pulls the foldseek fraction sequence identity (fident) from the Foldseek results files for each input protein.

    This will probably be replaced in the future by an all-v-all sequence identity comparison using FAMSA, WITCH, or other approach.
    """
    input:
        m8files=expand(
            output_dir / foldseekresults_dir / "{{protid}}" / "alis_{db}.m8", db=FS_DATABASES
        ),
    output:
        fident_features=output_dir / clusteringresults_dir / "{protid}_fident_features.tsv",
    params:
        protid="{protid}",
    benchmark:
        output_dir / benchmarks_dir / "{protid}.aggregate_foldseek_fident.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/aggregate_foldseek_fraction_seq_identity.py -i {input.m8files} -o {output.fident_features} -p {params.protid}
        """


rule aggregate_lists:
    """
    Take all Uniprot ID lists and make them one big ID list, removing duplicates.
    """
    input:
        expand(output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt", protid=PROTID),
        expand(output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt", protid=PROTID),
    output:
        jointlist=output_dir / clusteringresults_dir / "jointhits.txt",
    benchmark:
        output_dir / benchmarks_dir / "aggregate_lists.txt"
    shell:
        """
        python ProteinCartography/aggregate_lists.py -i {input} -o {output.jointlist}
        """


rule fetch_uniprot_metadata:
    """
    Use the output.jointlist file to query Uniprot and download all metadata as a big ol' TSV.
    """
    input:
        output_dir / clusteringresults_dir / "jointhits.txt",
    output:
        output_dir / clusteringresults_dir / "uniprot_features.tsv",
    params:
        additional_fields=UNIPROT_ADDITIONAL_FIELDS,
    benchmark:
        output_dir / benchmarks_dir / "get_uniprot_metadata.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py -i {input} -o {output} -a {params.additional_fields}
        """


rule filter_uniprot_hits:
    """
    Use the metadata features from Uniprot to filter hits based on sequence status, fragment, and size.
    """
    input:
        output_dir / clusteringresults_dir / "uniprot_features.tsv",
    output:
        output_dir / clusteringresults_dir / "alphafold_querylist.txt",
    params:
        min_length=MIN_LENGTH,
        max_length=MAX_LENGTH,
    benchmark:
        output_dir / benchmarks_dir / "filter_uniprot_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/filter_uniprot_hits.py -i {input} -o {output} -m {params.min_length} -M {params.max_length} --excluded-protids {PROTID}
        """


checkpoint download_pdbs:
    """
    Download all PDB files from AlphaFold
    """
    input:
        jointlist=output_dir / clusteringresults_dir / "alphafold_querylist.txt",
    output:
        output_dir=directory(pdb_download_dir),
    params:
        max_structures=MAX_STRUCTURES,
    benchmark:
        output_dir / benchmarks_dir / "download_pdbs.txt"
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py -i {input.jointlist} -o {output.output_dir} -M {params.max_structures}
        """


def checkpoint_download_pdbs(wildcards):
    """
    Returns the paths to all PDB files
    (both those downloaded in `download_pdbs` and those corresponding to the input proteins)
    """
    # the directory containing the PDB files
    pdb_dirpath = checkpoints.download_pdbs.get(**wildcards).output.output_dir

    # the paths to the PDB files downloaded in `download_pdbs`
    pdb_filepaths = sorted(Path(pdb_dirpath).glob("*.pdb"))

    # append the paths to the PDB files corresponding to the input proteins
    # note: this triggers the `copy_pdb` rule to copy the input PDB files from `input_dir`
    pdb_filepaths += expand(pdb_download_dir / "{protid}.pdb", protid=PROTID)
    return pdb_filepaths


rule assess_pdbs:
    """
    Calculates the quality of all PDBs
    (those downloaded from AlphaFold and those corresponding to the input proteins).
    """
    input:
        checkpoint_download_pdbs,
    output:
        pdb_features=output_dir / clusteringresults_dir / "pdb_features.tsv",
    params:
        pdb_download_dir=pdb_download_dir,
    benchmark:
        output_dir / benchmarks_dir / "assess_pdbs.txt"
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py -i {params.pdb_download_dir} -o {output.pdb_features}
        """


rule foldseek_clustering:
    """
    Runs foldseek all-v-all TM-score comparison and foldseek clustering.
    """
    input:
        checkpoint_download_pdbs,
    output:
        allvall_pivot=output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
        struclusters_features=output_dir / clusteringresults_dir / "struclusters_features.tsv",
    params:
        querydir=pdb_download_dir,
        resultsdir=output_dir / clusteringresults_dir,
    conda:
        "envs/foldseek.yml"
    resources:
        mem_mb=32 * 1000,
    threads: 16
    benchmark:
        output_dir / benchmarks_dir / "foldseek_clustering.txt"
    shell:
        """
        python ProteinCartography/foldseek_clustering.py -q {params.querydir} -r {params.resultsdir}
        """


rule dim_reduction:
    """
    Perform dimensionality reduction, saving as an embedding matrix and a TSV
    Write a set of functions to return Dataframes for interactive compute
    Write helper functions to save the dataframes only called by main()
    """
    input:
        output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
    output:
        output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted_{modes}.tsv",
    params:
        modes="{modes}",
    conda:
        "envs/analysis.yml"
    benchmark:
        output_dir / benchmarks_dir / "{modes}.dim_reduction.txt"
    shell:
        """
        python ProteinCartography/dim_reduction.py -i {input} -m {params.modes}
        """


rule leiden_clustering:
    """
    Performs Leiden clustering on the data using scanpy's implementation.
    """
    input:
        output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
    output:
        output_dir / clusteringresults_dir / "leiden_features.tsv",
    conda:
        "envs/analysis.yml"
    benchmark:
        output_dir / benchmarks_dir / "leiden_clustering.txt"
    shell:
        """
        python ProteinCartography/leiden_clustering.py -i {input} -o {output}
        """


rule input_distances:
    """
    Extracts the distances from input proteins to other proteins in the dataset.
    Adds them as options for the visualization plot.
    """
    input:
        output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
    output:
        output_dir / clusteringresults_dir / "{protid}_distance_features.tsv",
    params:
        protid="{protid}",
    benchmark:
        output_dir / benchmarks_dir / "{protid}.input_distances.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/extract_input_distances.py -i {input} -o {output} -p {params.protid}
        """


rule calculate_concordance:
    """
    Currently, this subtracts the fraction sequence identity from the TM-score to get a measure of whether something is more similar in sequence or structure.

    We're working on developing some kind of test statistic that evaluates the significance of this difference from some expectation.
    """
    input:
        tmscore_file=output_dir / clusteringresults_dir / "{protid}_distance_features.tsv",
        fident_file=output_dir / clusteringresults_dir / "{protid}_fident_features.tsv",
    params:
        protid="{protid}",
    output:
        output_dir / clusteringresults_dir / "{protid}_concordance_features.tsv",
    benchmark:
        output_dir / benchmarks_dir / "{protid}.calculate_concordance.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/calculate_concordance.py -t {input.tmscore_file} -f {input.fident_file} -o {output} -p {params.protid}
        """


rule get_source:
    """
    Checks the blasthits and foldseekhits files to determine the source of each protein.
    Adds this info to the visualization plot.
    """
    input:
        pivoted=output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
        hitfiles=expand(
            output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt", protid=PROTID
        )
        + expand(output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt", protid=PROTID),
    output:
        pivot=output_dir / clusteringresults_dir / "source_features.tsv",
    params:
        protid=PROTID,
    benchmark:
        output_dir / benchmarks_dir / "get_source.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/get_source.py -i {input.pivoted} -f {input.hitfiles} -o {output.pivot} -k {params.protid}
        """


rule aggregate_features:
    """
    Aggregate all TSV features provided by user in some specific directory, making one big TSV
    """
    input:
        output_dir / clusteringresults_dir / "struclusters_features.tsv",
        output_dir / clusteringresults_dir / "uniprot_features.tsv",
        output_dir / clusteringresults_dir / "leiden_features.tsv",
        expand(output_dir / clusteringresults_dir / "{protid}_distance_features.tsv", protid=PROTID),
        expand(output_dir / clusteringresults_dir / "{protid}_fident_features.tsv", protid=PROTID),
        expand(
            output_dir / clusteringresults_dir / "{protid}_concordance_features.tsv", protid=PROTID
        ),
        output_dir / clusteringresults_dir / "source_features.tsv",
        output_dir / clusteringresults_dir / "pdb_features.tsv",
    output:
        output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features.tsv"),
    params:
        override=OVERRIDE_FILE,
    benchmark:
        output_dir / benchmarks_dir / "aggregate_features.txt"
    conda:
        "envs/pandas.yml"
    shell:
        """
        python ProteinCartography/aggregate_features.py -i {input} -o {output} -v {params.override}
        """


rule plot_interactive:
    """
    Generate interactive scatter plot HTML programmatically based on user-input parameters
    Takes the TSV from rule aggregate_features and select default columns
    User should be able to call this module and pass their own functions to parse particular TSV columns
    Should have means to set a palette for each individual plot type, maybe as JSON?
    """
    input:
        dimensions=output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted_{modes}.tsv",
        features=output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features.tsv"),
    output:
        output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features_{modes}.html"),
    params:
        modes="{modes}",
        protid=expand("{protid}", protid=PROTID),
        taxon_focus=TAXON_FOCUS,
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "{modes}.plot_interactive.txt"
    shell:
        """
        python ProteinCartography/plot_interactive.py -d {input.dimensions} -f {input.features} -o {output} -t {params.modes} -k {params.protid} -x {params.taxon_focus}
        """


rule plot_similarity_leiden:
    """
    Plots a similarity score matrix for Leiden clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.
    """
    input:
        matrix=output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
        features=output_dir / clusteringresults_dir / "leiden_features.tsv",
    output:
        tsv=output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.tsv"),
        html=output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.html"),
    params:
        column="LeidenCluster",
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_similarity_leiden.txt"
    shell:
        """
        python ProteinCartography/cluster_similarity.py -m {input.matrix} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        """


rule plot_similarity_strucluster:
    """
    Plots a similarity score matrix for Foldseek's structural clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.
    """
    input:
        matrix=output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted.tsv",
        features=output_dir / clusteringresults_dir / "struclusters_features.tsv",
    output:
        tsv=output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.tsv"),
        html=output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.html"),
    params:
        column="StruCluster",
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_similarity_strucluster.txt"
    shell:
        """
        python ProteinCartography/cluster_similarity.py -m {input.matrix} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        """


rule plot_semantic_analysis:
    """
    Plots a semantic analysis chart for groups within the data.
    """
    input:
        features_file=output_dir
        / clusteringresults_dir
        / (analysis_name + "_aggregated_features.tsv"),
    output:
        pdf=output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.pdf"),
        interactive=output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.html"),
    params:
        agg_column="LeidenCluster",
        annot_column="'Protein names'",
        analysis_name=analysis_name,
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_semantic_analysis.txt"
    shell:
        """
        python ProteinCartography/semantic_analysis.py -f {input.features_file} -c {params.agg_column} -n {params.annot_column} -o {output.pdf} -i {output.interactive} -a {params.analysis_name}
        """


rule plot_cluster_distributions:
    """
    Plots distributions of key values per cluster for each input protein.
    """
    input:
        features_file=output_dir
        / clusteringresults_dir
        / (analysis_name + "_aggregated_features.tsv"),
    output:
        output_dir / clusteringresults_dir / (analysis_name + "_{protid}_distribution_analysis.svg"),
    params:
        protid="{protid}",
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_cluster_distributions_{protid}.txt"
    shell:
        """
        python ProteinCartography/plot_cluster_distributions.py -i {input.features_file} -o {output} -k {params.protid}
        """
