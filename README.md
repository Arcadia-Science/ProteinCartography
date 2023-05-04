# gene-family-cartography
embedding proteins into feature spaces for computational exploration

## Purpose
The relationship between protein sequence, structure, and function has only been thoroughly investigated for a handful of gene families.  
This repo takes an agnostic approach to characterizing groups of similar proteins using feature embedding spaces.  

Our pipeline starts with user-provided protein(s) of interest and searches the available sequence and structure databases for matches.  
Using the full list of matches, we can build a "map" of all the similar proteins and look for clusters of proteins with similar features.  
Overlaying a variety of different parameters such as taxonomy, sequence divergence, and other features onto these spaces allows us to explore the features that drive differences between clusters.

## Directory Structure
- [Snakefile](Snakefile): the Snakemake pipeline that orchestrates this repo's functions.
- [cartography.yml](cartography.yml): the conda environment for this repo.
- [utils/](utils/): helper scripts that are called by the Snakemake pipeline, along with a few other bits and bobs
- [notebooks/](notebooks/): Jupyter notebook examples of manually running this data. These are messy and will be tidied up in the final repo. Currently a development playground - don't look too deeply here yet.
- [examples/](examples/): example output files.

## Pipeline Overview
This repo builds a Snakemake pipeline that takes input `.fasta` and `.pdb` files.
The rulegraph for this pipeline is as follows:  
![rulegraph](rulegraph.png)

The steps of the pipeline have the following functionality:

### Protein Folding
0. Fold input FASTA files using ESMFold.
    - The workflow starts with input FASTA files, one entry per file, with a unique ID.
    - If a matching PDB file is not provided, the pipeline will use the ESMFold API to generate a PDB.

### Protein Search
1. Search the Alphafold databases using queries to the FoldSeek webserver API for each provided `.pdb` file.  
    - Currently, we're limited to a maximum of 1000 hits per protein.  
    - We should look for ways to increase this number of hits, potentially by hosting the full AlphaFold database ourselves and building a simple API to query it.
    - This fails more often than we'd like because Foldseek rate limits our requests after a certain threshold. The actual threshold is unknown.  
2. Search the non-redundant GenBank/RefSeq database using blastp for each provided `.fasta` file.  
    - Takes the resulting output hits and maps each GenBank/RefSeq hit to a Uniprot ID using [**`bioservices UniProt`**](https://bioservices.readthedocs.io/en/latest/references.html#module-bioservices.uniprot)
    - This also fails more often than we'd like if multiple queries are made in a short time period.

### Download Data
3. Aggregate the list of Foldseek and BLAST hits into a single list of Uniprot IDs.  

4. Download a `.pdb` file from each protein from AlphaFold.  
    - This part is a bit slow, as it's currently limited to the number of cores provided to Snakemake.

---
== TODO ==

5. Download annotation and feature information for each protein from Uniprot.  
    - This pretty easy using **`bioservices UniProt`**.  

### Clustering

6. Cluster all protein .pdb files using FoldSeek.  
7. Generate a similarity matrix.

### Data Aggregation

8. Take features from Uniprot and other manually-input sources and aggregate them into one large .tsv file.
    - Will need to build this in a modular fashion.

### Plotting

9. Build an explorable HTML visualization using `Plotly` based on the aggregated features.  
    - An example can be found [here](examples/scatter.html)
    - Each point has hover-over information
    - Default parameters would include:
        - Foldseek Structural cluster
        - Leiden cluster
        - Annotation score
        - Protein length
    - Power users would be able to customize the visualization using a variety of built-in parameters.

== FUTURE ==
> Processes below haven't been explored yet but are of interest.

- Generate distance matrices using sequence similarity instead of structure.
- Allow for passing of arbitrary TSV data types for building visualizations.
- Automated aggregation of input/output files to share using `biofile`?
