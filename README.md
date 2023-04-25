# gene-family-cartography
a pipeline to embed gene families into distance spaces for exploration

## Overview
This repo builds a Snakemake pipeline that takes input `.fasta` and `.pdb` files and does the following operations:

### Protein Search
1. Search the Alphafold databases using queries to the FoldSeek webserver API for each provided `.pdb` file.  
    - Currently, we're limited to a maximum of 1000 hits per protein.  
    - We should look for ways to increase this number of hits, potentially by hosting the full AlphaFold database ourselves and building a simple API to query it.  
2. Search the non-redundant GenBank/RefSeq database using blastp for each provided `.fasta` file.  
    - Takes the resulting output hits and maps each GenBank/RefSeq hit to a Uniprot ID using [**`bioservices UniProt`**](https://bioservices.readthedocs.io/en/latest/references.html#module-bioservices.uniprot)  

### Download Data
3. Aggregate the list of Foldseek and BLAST hits into a single list of Uniprot IDs.  

---
== TODO ==
> Processes below are tractable but haven't yet been built into the Snakemake pipeline.  

4. Download a `.pdb` file from each protein from AlphaFold.  
    - Need to figure out a way to do this that's highly parallelized to be as fast as possible.  

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
