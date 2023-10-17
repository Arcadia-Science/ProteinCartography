import os
from pathlib import Path

###########################################
## Parse config information
###########################################

# Default pipeline configuration parameters are in this file
# If you create a new yml file and use the --configfile flag,
# options in that new file overwrite the defaults.
configfile: './config.yml'

# Set the input directory
#### In the future, also accept Uniprot accession numbers, which will be auto-queried and downloaded
input_dir = Path(config["input_dir"])

# Set the prefix of the output file of the analysis
analysis_name = config["analysis_name"]

# put most things into the output directory
output_dir = Path(config["output_dir"])

# Check for an override file, setting a variable if it exists
if "override_file" in config:
    OVERRIDE_FILE = input_dir / config["override_file"]

    # If it isn't a real file, ignore it
    if not os.path.exists(OVERRIDE_FILE):
        OVERRIDE_FILE = ''
else:
    OVERRIDE_FILE = ''

if "taxon_focus" in config:
    TAXON_FOCUS = config["taxon_focus"]
else:
    TAXON_FOCUS = 'euk'

UNIPROT_ADDITIONAL_FIELDS = config["uniprot_additional_fields"]

MAX_BLASTHITS = int(config["max_blasthits"])
MAX_STRUCTURES = int(config["max_structures"])

FS_DATABASES = config["foldseek_databases"]
MODES = config["plotting_modes"]

MIN_LENGTH = int(config["min_length"])
MAX_LENGTH = int(config["max_length"])

###########################################
## Setup directory structure
###########################################

# these directories fall within the output directory
blastresults_dir = Path('blastresults/')
foldseekresults_dir = Path('foldseekresults/')
foldseekclustering_dir = Path('foldseekclustering/')
clusteringresults_dir = Path('clusteringresults/')
benchmarks_dir = Path('benchmarks/')
downloading_dir = Path('downloading/')

# gets the protein ID based on FASTA file name
# flexibly checks if fasta file is correct suffix
FASTA_FORMATS = ['.fa', '.fasta', '.fna']
PROTID = []
for file in os.listdir(input_dir):
    if any(file.lower().endswith(suffix) for suffix in FASTA_FORMATS):
        file_id = os.path.splitext(file)[0]
        PROTID.append(file_id)

BLAST_DEFAULTS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sacc', 'saccver', 'sgi', 'staxids', 'scomnames']
BLAST_DEFAULT_STRING = '"' + ' '.join(['6'] + BLAST_DEFAULTS) + '"'

######################################

rule all:
    input:
        expand(output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features_{modes}.html"), modes = MODES),
        output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.html"),
        output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.html"),
        output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.pdf"),
        output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.html")

###########################################
## make .pdb files using esmfold API query
###########################################

rule make_pdb:
    '''
    Use the ESMFold API query to generate a pdb from a fasta file.
    This rule is only touched if a .pdb file doesn't already exist.
    '''
    input:
        cds = input_dir / "{protid}.fasta"
    output:
        pdb = input_dir / "{protid}.pdb"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.make_pdb.txt"
    conda:
        "envs/api.yml"
    shell:
        '''
        python ProteinCartography/esmfold_apiquery.py -i {input.cds}
        '''

# Prioritize copying an existing PDB to the folder if the ID is also a blast/foldseek hit.
ruleorder: copy_pdb > download_pdbs

rule copy_pdb:
    '''
    Copies existing or generated PDBs to the Foldseek clustering folder.
    '''
    input: input_dir / "{protid}.pdb"
    output: output_dir / foldseekclustering_dir / "{protid}.pdb"
    shell:
        '''
        cp {input} {output}
        '''

###########################################
## perform blastp to full database using nr
###########################################

rule run_blast:
    '''
    Using files located in the input directory, perform BLAST using the web API.

    Large proteins will cause remote BLAST to fail; you can still perform a manual BLAST search to get around this.
    '''
    input:
        cds = input_dir / "{protid}.fasta"
    output:
        blastresults = output_dir / blastresults_dir / "{protid}.blastresults.tsv",
        refseqhits = output_dir / blastresults_dir / "{protid}.blasthits.refseq.txt"
    params:
        max_blasthits = MAX_BLASTHITS,
        blast_string = BLAST_DEFAULT_STRING
    benchmark:
        output_dir / benchmarks_dir / "{protid}.run_blast.txt"
    conda:
        "envs/run_blast.yml"
    shell:
        '''
        blastp -db nr -query {input.cds} -out {output.blastresults} -remote -max_target_seqs {params.max_blasthits} -outfmt {params.blast_string}
        python ProteinCartography/extract_blasthits.py -i {output.blastresults} -o {output.refseqhits} -B {params.blast_string}
        '''

rule map_refseqids:
    '''
    Using List of RefSeq IDs, query the Uniprot ID mapping tool using bioservices UniProt.
    Returns a list of UniProt IDs.
    '''
    input:
        refseqhits = output_dir / blastresults_dir / "{protid}.blasthits.refseq.txt"
    output:
        uniprothits = output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.map_refseqids.txt"
    conda:
        "envs/map_refseq_ids.yml"
    shell:
        '''
        python ProteinCartography/map_refseqids.py -i {input.refseqhits} -o {output.uniprothits}
        '''

######################################
## perform Foldseek using web API
######################################

# Note: this returns a limited number of hits; up to 1000 per database

rule run_foldseek:
    '''
    Runs Foldseek using a query to the web API using a custom Python script.
    The script accepts an input file ending in '.pdb' and returns an output file ending in '.tar.gz'.
    The script also accepts a `--mode` flag of either '3diaa' (default) or 'tmalign' and choice of databases.
    After running, untars the files and extracts hits.
    '''
    input:
        cds = input_dir / "{protid}.pdb"
    output:
        targz = output_dir / foldseekresults_dir / "{protid}.fsresults.tar.gz",
        unpacked = directory(output_dir / foldseekresults_dir / "{protid}"),
        m8files = expand(output_dir / foldseekresults_dir / "{{protid}}" / "alis_{db}.m8", db = FS_DATABASES),
        foldseekhits = output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt"
    params:
        fs_databases = expand("{fs_databases}", fs_databases = FS_DATABASES)
    conda:
        "envs/foldseek.yml"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.run_foldseek.txt"
    shell:
        '''
        python ProteinCartography/foldseek_apiquery.py -i {input.cds} -o {output.targz} -d {params.fs_databases}
        tar -xvf {output.targz} -C {output.unpacked}
        python ProteinCartography/extract_foldseekhits.py -i {output.m8files} -o {output.foldseekhits}
        '''

rule aggregate_foldseek_fraction_seq_identity:
    '''
    Pulls the foldseek fraction sequence identity (fident) from the Foldseek results files for each input protein.

    This will probably be replaced in the future by an all-v-all sequence identity comparison using FAMSA, WITCH, or other approach.
    '''
    input:
        m8files = expand(output_dir / foldseekresults_dir / "{{protid}}" / "alis_{db}.m8", db = FS_DATABASES)
    output:
        fident_features = output_dir / clusteringresults_dir / "{protid}_fident_features.tsv"
    params:
        protid = "{protid}"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.aggregate_foldseek_fident.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/aggregate_foldseek_fraction_seq_identity.py -i {input.m8files} -o {output.fident_features} -p {params.protid}
        '''

#####################################################################
## aggregate all hits, download metadata, download structure files
#####################################################################

rule aggregate_lists:
    '''
    Take all Uniprot ID lists and make them one big ID list, removing duplicates.
    '''
    input:
        expand(output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt", protid = PROTID),
        expand(output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt", protid = PROTID)
    output:
        jointlist = output_dir / clusteringresults_dir / "jointhits.txt"
    benchmark:
        output_dir / benchmarks_dir / "aggregate_lists.txt"
    shell:
        '''
        python ProteinCartography/aggregate_lists.py -i {input} -o {output.jointlist}
        '''

rule fetch_uniprot_metadata:
    '''
    Use the output.jointlist file to query Uniprot and download all metadata as a big ol' TSV.
    '''
    input: output_dir / clusteringresults_dir / "jointhits.txt"
    output: output_dir / clusteringresults_dir / "uniprot_features.tsv"
    params:
        additional_fields = UNIPROT_ADDITIONAL_FIELDS
    benchmark:
        output_dir / benchmarks_dir / "get_uniprot_metadata.txt"
    conda:
        "envs/fetch_uniprot_metadata.yml"
    shell:
        '''
        python ProteinCartography/fetch_uniprot_metadata.py -i {input} -o {output} -a {params.additional_fields}
        '''

rule filter_uniprot_hits:
    '''
    Use the metadata features from Uniprot to filter hits based on sequence status, fragment, and size.
    '''
    input: output_dir / clusteringresults_dir / "uniprot_features.tsv"
    output: output_dir / clusteringresults_dir / "alphafold_querylist.txt"
    params:
        min_length = MIN_LENGTH,
        max_length = MAX_LENGTH
    benchmark:
        output_dir / benchmarks_dir / "filter_uniprot_hits.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/filter_uniprot_hits.py -i {input} -o {output} -m {params.min_length} -M {params.max_length}
        '''

checkpoint create_alphafold_wildcard:
    '''
    Create dummy files to make Snakemake detect a wildcard.
    '''
    input:
        jointlist = output_dir / clusteringresults_dir / "alphafold_querylist.txt"
    output: directory(os.path.join(output_dir, "alphafold_dummy/"))
    params:
        max_structures = MAX_STRUCTURES
    shell:
        '''
        python ProteinCartography/make_dummies.py -i {input.jointlist} -o {output} -M {params.max_structures}
        '''

rule download_pdbs:
    '''
    Use a checkpoint to parse all of the items in the output.jointlist file from aggregate lists and download all the PDBs.
    '''
    input:
        output_dir / "alphafold_dummy/{acc}.txt"
    output:
        output_dir / foldseekclustering_dir / "{acc}.pdb"
    params:
        outdir = output_dir / foldseekclustering_dir
    benchmark:
        output_dir / benchmarks_dir / downloading_dir / "{acc}.download_pdbs.txt"
    conda:
        "envs/api.yml"
    resources:
        mem_mb=256
    threads: 1
    localrule: True
    shell:
        '''
        python ProteinCartography/fetch_accession.py -a {wildcards.acc} -o {params.outdir} -f pdb
        '''

def checkpoint_create_alphafold_wildcard(wildcards):
    # expand checkpoint to get acc values
    checkpoint_output = checkpoints.create_alphafold_wildcard.get(**wildcards).output[0]

    # trawls the checkpoint_output file for .txt files
    # and generates expected .pdb file names for foldseekclustering_dir
    file_names = expand(output_dir / foldseekclustering_dir / "{acc}.pdb",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}.txt")).acc) + \
                 expand(output_dir / foldseekclustering_dir / "{protid}.pdb", protid = PROTID)
    return file_names

rule assess_pdbs:
    '''
    Calculates the quality of all PDBs downloaded from AlphaFold.
    '''
    input: checkpoint_create_alphafold_wildcard
    output:
        pdb_paths = output_dir / clusteringresults_dir / 'pdbpaths.txt',
        pdb_features = output_dir / clusteringresults_dir / 'pdb_features.tsv'
    params:
        inputdir = input_dir,
        clusteringdir = output_dir / foldseekclustering_dir
    benchmark:
        output_dir / benchmarks_dir / "assess_pdbs.txt"
    conda:
        "envs/plotting.yml"
    shell:
        '''
        python ProteinCartography/prep_pdbpaths.py -d {params.clusteringdir} {params.inputdir} -o {output.pdb_paths}
        python ProteinCartography/assess_pdbs.py -t {output.pdb_paths} -o {output.pdb_features}
        '''

#####################################################################
## clustering and dimensionality reduction
#####################################################################

rule foldseek_clustering:
    '''
    Runs foldseek all-v-all TM-score comparison and foldseek clustering.
    '''
    input: checkpoint_create_alphafold_wildcard
    output:
        allvall_pivot = output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv',
        struclusters_features = output_dir / clusteringresults_dir / 'struclusters_features.tsv'
    params:
        querydir = output_dir / foldseekclustering_dir,
        resultsdir = output_dir / clusteringresults_dir
    conda:
        "envs/foldseek.yml"
    resources:
        mem_mb=32 * 1000
    threads: 16
    benchmark:
        output_dir / benchmarks_dir / "foldseek_clustering.txt"
    shell:
        '''
        python ProteinCartography/foldseek_clustering.py -q {params.querydir} -r {params.resultsdir}
        '''

rule dim_reduction:
    '''
    Perform dimensionality reduction, saving as an embedding matrix and a TSV
    Write a set of functions to return Dataframes for interactive compute
    Write helper functions to save the dataframes only called by main()
    '''
    input: output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv'
    output: output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted_{modes}.tsv'
    params:
        modes = '{modes}'
    conda:
        "envs/dim_reduction.yml"
    benchmark:
        output_dir / benchmarks_dir / "{modes}.dim_reduction.txt"
    shell:
        '''
        python ProteinCartography/dim_reduction.py -i {input} -m {params.modes}
        '''

rule leiden_clustering:
    '''
    Performs Leiden clustering on the data using scanpy's implementation.
    '''
    input: output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv'
    output: output_dir / clusteringresults_dir / 'leiden_features.tsv'
    conda:
        "envs/leiden_clustering.yml"
    benchmark:
        output_dir / benchmarks_dir / "leiden_clustering.txt"
    shell:
        '''
        python ProteinCartography/leiden_clustering.py -i {input} -o {output}
        '''

rule input_distances:
    '''
    Extracts the distances from input proteins to other proteins in the dataset.
    Adds them as options for the visualization plot.
    '''
    input: output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv'
    output: output_dir / clusteringresults_dir / '{protid}_distance_features.tsv'
    params:
        protid = "{protid}"
    benchmark:
        output_dir / benchmarks_dir / "{protid}.input_distances.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/extract_input_distances.py -i {input} -o {output} -p {params.protid}
        '''

rule calculate_concordance:
    '''
    Currently, this subtracts the fraction sequence identity from the TM-score to get a measure of whether something is more similar in sequence or structure.

    We're working on developing some kind of test statistic that evaluates the significance of this difference from some expectation.
    '''
    input:
        tmscore_file = output_dir / clusteringresults_dir / '{protid}_distance_features.tsv',
        fident_file = output_dir / clusteringresults_dir / '{protid}_fident_features.tsv'
    params:
        protid = "{protid}"
    output: output_dir / clusteringresults_dir / '{protid}_concordance_features.tsv'
    benchmark:
        output_dir / benchmarks_dir / "{protid}.calculate_concordance.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/calculate_concordance.py -t {input.tmscore_file} -f {input.fident_file} -o {output} -p {params.protid}
        '''

rule get_source:
    '''
    Checks the blasthits and foldseekhits files to determine the source of each protein.
    Adds this info to the visualization plot.
    '''
    input:
        pivoted = output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv',
        hitfiles = expand(output_dir / foldseekresults_dir / "{protid}.foldseekhits.txt", protid = PROTID) + expand(output_dir / blastresults_dir / "{protid}.blasthits.uniprot.txt", protid = PROTID)
    output:
        pivot = output_dir / clusteringresults_dir / 'source_features.tsv'
    params:
        protid = PROTID
    benchmark:
        output_dir / benchmarks_dir / "get_source.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/get_source.py -i {input.pivoted} -f {input.hitfiles} -o {output.pivot} -k {params.protid}
        '''

#####################################################################
## aggregate features into a big TSV and make a nice plot
#####################################################################

rule aggregate_features:
    '''
    Aggregate all TSV features provided by user in some specific directory, making one big TSV
    '''
    input:
        output_dir / clusteringresults_dir / "struclusters_features.tsv",
        output_dir / clusteringresults_dir / "uniprot_features.tsv",
        output_dir / clusteringresults_dir / "leiden_features.tsv",
        expand(output_dir / clusteringresults_dir / "{protid}_distance_features.tsv", protid = PROTID),
        expand(output_dir / clusteringresults_dir / "{protid}_fident_features.tsv", protid = PROTID),
        expand(output_dir / clusteringresults_dir / "{protid}_concordance_features.tsv", protid = PROTID),
        output_dir / clusteringresults_dir / "source_features.tsv",
        output_dir / clusteringresults_dir / "pdb_features.tsv"
    output: output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features.tsv")
    params:
        override = OVERRIDE_FILE
    benchmark:
        output_dir / benchmarks_dir / "aggregate_features.txt"
    conda:
        "envs/pandas.yml"
    shell:
        '''
        python ProteinCartography/aggregate_features.py -i {input} -o {output} -v {params.override}
        '''

rule plot_interactive:
    '''
    Generate interactive scatter plot HTML programmatically based on user-input parameters
    Takes the TSV from rule aggregate_features and select default columns
    User should be able to call this module and pass their own functions to parse particular TSV columns
    Should have means to set a palette for each individual plot type, maybe as JSON?
    '''
    input:
        dimensions = output_dir / clusteringresults_dir / "all_by_all_tmscore_pivoted_{modes}.tsv",
        features = output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features.tsv")
    output:
        output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features_{modes}.html")
    params:
        modes = "{modes}",
        protid = expand("{protid}", protid = PROTID),
        taxon_focus = TAXON_FOCUS
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "{modes}.plot_interactive.txt"
    shell:
        '''
        python ProteinCartography/plot_interactive.py -d {input.dimensions} -f {input.features} -o {output} -t {params.modes} -k {params.protid} -x {params.taxon_focus}
        '''

rule plot_similarity_leiden:
    '''
    Plots a similarity score matrix for Leiden clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.
    '''
    input:
        matrix = output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv',
        features = output_dir / clusteringresults_dir / "leiden_features.tsv",
    output:
        tsv = output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.tsv"),
        html = output_dir / clusteringresults_dir / (analysis_name + "_leiden_similarity.html")
    params:
        column = 'LeidenCluster'
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_similarity_leiden.txt"
    shell:
        '''
        python ProteinCartography/cluster_similarity.py -m {input.matrix} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        '''

rule plot_similarity_strucluster:
    '''
    Plots a similarity score matrix for Foldseek's structural clusters.
    For each cluster, calculates the mean TM-score of all structures in that cluster versus all other clusters.
    The diagonal of the plot shows how similar proteins are within a given cluster.
    The other cells show how similar other clusters are to each other.
    '''
    input:
        matrix = output_dir / clusteringresults_dir / 'all_by_all_tmscore_pivoted.tsv',
        features = output_dir / clusteringresults_dir / "struclusters_features.tsv",
    output:
        tsv = output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.tsv"),
        html = output_dir / clusteringresults_dir / (analysis_name + "_strucluster_similarity.html")
    params:
        column = 'StruCluster'
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_similarity_strucluster.txt"
    shell:
        '''
        python ProteinCartography/cluster_similarity.py -m {input.matrix} -f {input.features} -c {params.column} -T {output.tsv} -H {output.html}
        '''

rule plot_semantic_analysis:
    '''
    Plots a semantic analysis chart for groups within the data.
    '''
    input:
        features_file = output_dir / clusteringresults_dir / (analysis_name + "_aggregated_features.tsv")
    output:
        pdf = output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.pdf"),
        interactive = output_dir / clusteringresults_dir / (analysis_name + "_semantic_analysis.html")
    params:
        agg_column = 'LeidenCluster',
        annot_column = "'Protein names'",
        analysis_name = analysis_name
    conda:
        "envs/plotting.yml"
    benchmark:
        output_dir / benchmarks_dir / "plot_semantic_analysis.txt"
    shell:
        '''
        python ProteinCartography/semantic_analysis.py -f {input.features_file} -c {params.agg_column} -n {params.annot_column} -o {output.pdf} -i {output.interactive} -a {params.analysis_name}
        '''
