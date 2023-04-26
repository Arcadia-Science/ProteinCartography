import os

# Get protein IDs from fasta files in input directory
# In the future, have this specified by command line argument or config file
# In the future, also accept Uniprot accession numbers, which will be auto-queried and downloaded
inputdir = 'input/'
PROTID = [os.path.basename(i).split('.fasta')[0] for i in os.listdir(inputdir) if '.fasta' in i]

######################################

rule all:
    input:
        expand("output/blastresults/{protid}.blasthits.uniprot.txt", protid = PROTID),
        expand("output/foldseekresults/{protid}.foldseekhits.txt", protid = PROTID),
        "output/foldseekclustering/alphafold_querylist.txt"

###########################################
## make .pdb files using gget alphafold
###########################################

# Ignore this for now, I'm having hardware problems getting gget to work
# It seems to be an Apple M1-specific problem

rule make_pdb:
    '''
    Use gget alphafold to generate a pdb from a fasta file.
    '''
    input:
        cds = "input/{protid}.fasta"
    output:
        pdb = "input/{protid}.pdb"
    shell:
        '''
        gget alphafold {input.cds} > {output.pdb}
        '''

###########################################
## perform blastp to full database using nr
###########################################

rule run_blast:
    '''
    Using files located in the "inputs/" folder, perform BLAST using the web API.
    '''
    input:
        cds = "input/{protid}.fasta"
    output:
        results = "output/blastresults/{protid}.blastresults.tsv"
    shell:
        '''
        blastp -db nr -query {input.cds} -out {output.results} -remote -max_target_seqs 50000 -outfmt 6
        '''

rule extract_blasthits:
    '''
    Using blast results files, generate lists of RefSeq ids for ID mapping.
    '''
    input:
        blastresults = "output/blastresults/{protid}.blastresults.tsv"
    output:
        refseqhits = "output/blastresults/{protid}.blasthits.refseq.txt"
    shell:
        '''
        python utils/extract_blasthits.py -i {input.blastresults} -o {output.refseqhits}
        '''

rule map_refseqids:
    '''
    Using List of RefSeq IDs, query the Uniprot ID mapping tool using bioservices UniProt.
    Returns a list of UniProt IDs.
    '''
    input:
        refseqhits = "output/blastresults/{protid}.blasthits.refseq.txt"
    output:
        uniprothits = "output/blastresults/{protid}.blasthits.uniprot.txt"
    shell:
        '''
        python utils/map_refseqids.py -i {input.refseqhits} -o {output.uniprothits}
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
    '''
    input:
        cds = "input/{protid}.pdb"
    output:
        results = "output/foldseekresults/{protid}.fsresults.tar.gz"
    shell:
        '''
        python utils/foldseek_apiquery.py -i {input.cds} -o {output.results}
        '''

rule unpack_fsresults:
    '''
    Untars and unzips files into a directory for each protid.
    '''
    input:
        targz = "output/foldseekresults/{protid}.fsresults.tar.gz"
    output:
        unpacked = directory("output/foldseekresults/{protid}/")
    shell:
        '''
        mkdir {output.unpacked}
        tar -xvf {input.targz} -C {output.unpacked}
        '''

rule extract_foldseekhits:
    '''
    Using Foldseek results directory, generate lists of Uniprot ids for mapping to Uniprot and pulling down PDB files.
    '''
    input:
        foldseekresults = "output/foldseekresults/{protid}/"
    output:
        foldseekhits = "output/foldseekresults/{protid}.foldseekhits.txt"
    shell:
        '''
        python utils/extract_foldseekhits.py -i {input.foldseekresults} -o {output.foldseekhits}
        '''

#####################################################################
## aggregate all hits, download metadata, download structure files
#####################################################################

rule aggregate_lists:
    '''
    Take all Uniprot ID lists and make them one big ID list, removing duplicates.
    '''
    input:
        foldseekresults_dir = "output/foldseekresults/",
        blastresults_dir = "output/blastresults/"
    params:
        foldseekresults_suffix = ".foldseekhits.txt",
        blastresults_suffix = ".blasthits.uniprot.txt"
    output:
        jointlist = "output/foldseekclustering/alphafold_querylist.txt"
    shell:
        '''
        python utils/aggregate_lists.py -i {input.foldseekresults_dir} {input.blastresults_dir} \
        -s {params.foldseekresults_suffix} {params.blastresults_suffix} \
        -o {output.jointlist}
        '''
        