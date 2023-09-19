# ProteinCartography
The ProteinCartography pipeline searches sequence and structure databases for matches to input proteins and builds maps of protein space for the purposes of discovery and exploration.

---
## Purpose
The relationship between protein sequence, structure, and function has only been thoroughly investigated for a handful of gene families. This repo takes an agnostic approach to characterizing groups of similar proteins using feature embedding spaces.  

Our pipeline starts with user-provided protein(s) of interest and searches the available sequence and structure databases for matches. Using the full list of matches, we can build a "map" of all the similar proteins and look for clusters of proteins with similar features. Overlaying a variety of different parameters such as taxonomy, sequence divergence, and other features onto these spaces allows us to explore the features that drive differences between clusters.

---
## Quickstart
1. Clone the GitHub repository.
    ```
    git clone https://github.com/Arcadia-Science/ProteinCartography.git
    ```
2. Install [`conda`](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and/or [`mamba`](https://github.com/mamba-org/mamba) if you don't already have them installed.
3. Create a conda environment containing the software needed to run the pipeline.  
    Run this code from within the ProteinCartography repo.
    ```
    conda env create -f envs/cartography_tidy.yml -n cartography_tidy
    conda activate cartography_tidy
    ```
4. Run the Snakemake pipeline using a demo protein (human ACTB, [P60709](https://www.uniprot.org/uniprotkb/P60709/entry)).
    Set `n` to be the number of cores you'd like to use for running the pipeline.
    ```
    snakemake --snakefile Snakefile --configfile demo/config_actin.yml --use-conda --cores n
    ```
5. Inspect results.
    In the `demo/output/clusteringresults/` directory, you should find the following files:  
    - `actin_aggregated_features.tsv`: metadata file containing protein feature hits
    - `actin_aggregated_features_pca_umap.html`: interactive UMAP scatter plot of results
    - `actin_aggregated_features_pca_tsne.html`: interactive t-SNE scatter plot of results
    - `actin_leiden_similarity.html`: mean cluster TM-score similarity heatmap
    - `actin_semantic_analysis.html` and `actin_semantic_analysis.pdf`: simple semantic analysis of clusters

---
## Modes
The pipeline supports two modes: **From-Query** and **From-Folder**.

### From-Query Mode
In this mode, the pipeline starts with a set of input proteins of interest in PDB and FASTA format and performs broad BLAST and Foldseek searches to identify hits. The pipeline aggregates all hits, downloads PDBs, and builds a map. The pipeline is implemented in [`Snakefile`](Snakefile).
![rulegraph](rulegraph.png)

#### Inputs
- Query protein files in FASTA and PDB format.
    - Each protein should have a unique protein identifier (`protid`).
    - The `protid` should be the prefix of the FASTA and PDB files (e.g. `P60709.fasta`, `P60709.pdb`).
    - Proteins provided as FASTAs <400aa in length will be folde automatically by the pipeline using the ESMFold API.
    - For proteins >400aa in length, we recommend using [ColabFold](https://github.com/sokrypton/ColabFold). 
- `config.yml` file with custom settings.
    - [`config.yml`](config.yml) is an example file that contains the defaults of the pipeline.
    - We recommend making a copy of this file and adjusting the settings.
    - A custom config file passed through the `--configfile` flag to Snakemake will overwrite the defaults.
    - **Key Parameters:**
        - `input`: directory containing input PDBs and FASTAs.
        - `output`: directory where all pipeline outputs are placed.
        - `analysis_name`: nickname for the analysis, appended to important output files.
        - See [`config.yml`](config.yml) for additional parameters.
- (Optional) `features_override.tsv` file.
    - This file can be used to pass protein metadata for custom proteins not found in UniProt.
    - These proteins will not automatically have metadata retrieved, so providing this information is useful for visualizing custom input proteins.
    - The columns of this file are described in [Feature file main columns](#feature-file-main-columns).
    - It can also be used to overwrite metadata values for specific proteins.

#### Walkthrough
0. Follow steps 1-3 in the [Quickstart](#quickstart) section above to set up and activate the conda environment.
1. Set up a [`config.yml`](config.yml) file specifying input and output directories and an analysis name.
2. Add query protein files in FASTA and PDB format to the input folder.
    - You can pull down proteins with a Uniprot accession number and Alphafold structure using the following command.
        Replace {accession} with your Uniprot accession (e.g. P24928).
        `python ProteinCartography/fetch_accession.py -a {accession} -o input -f fasta pdb`
        This saves a FASTA file from Uniprot and a PDB file from AlphaFold to the `input/` folder.
3. Run this command, replacing the `config.yml` with your config file and `8` with the number of cores you want to allocate to Snakemake.
```
snakemake --snakefile Snakefile --configfile config.yml --use-conda --cores 8
```

### From-Folder Mode
In this mode, the pipeline starts with a folder containing PDBs of interest and performs just the clustering and visualization steps of the pipeline, without performing any searches or downloads. The pipeline is implemented in [`Snakefile_ff`](Snakefile_ff).
![rulegraph_ff](rulegraph_ff.png)

#### Inputs
- Input protein structures in PDB format.
    - Each protein should have a unique protein identifier (`protid`) which matches the PDB file prefix, as described for [From-Query mode](#from-query-mode).
    - The pipeline does not yet support proteins with multiple chains.
- `config_ff.yml` file with custom settings.
    - [`config_ff.yml`](config_ff.yml) is an example file that contains the defaults of the pipeline.
    - We recommend making a copy of this file and adjusting the settings.
    - A custom config file passed through the `--configfile` flag to Snakemake will overwrite the defaults.
    - **Key Parameters:**
        - `input`: directory containing input PDBs and FASTAs.
        - `output`: directory where all pipeline outputs are placed.
        - `analysis_name`: nickname for the analysis, appended to important output files.
        - `features_file`: path to features file (described below).
        - (Optional) `keyids`: key `protid` values for proteins to highlight similar to input proteins in the From-Query mode.
        - See [`config_ff.yml`](config_ff.yml) for additional parameters.
- Features file with protein metadata.
    - Usually, we call this file `uniprot_features.tsv` but you can use any name.
    - The file should be placed into the `input` directory.
    - This file contains protein metadata for each protein used for visualization purposes.
    - The columns of this file are described in [Feature file main columns](#feature-file-main-columns).
    
#### Walkthrough
0. Follow steps 1-3 in the [Quickstart](#quickstart) section above to set up and activate the conda environment.
1. Set up a [`config_ff.yml`](config_ff.yml) file specifying input and output directories and an analysis name.
2. Add input protein structures in PDB format to your input folder.
    - The pipeline does not yet support PDBs with multiple chains.
3. Make a `uniprot_features.tsv` file.
    You can generate this one of two ways. If you have a list of Uniprot accessions, you can provide that as a .txt file and automatically pull down the correct annotations. Alternatively, you can manually generate the file.  
    **3.1 Generating the file using query_uniprot.py**
    - Create a `uniprot_ids.txt` file that contains a list of Uniprot accessions, one per line.
    - Run the following code from the base of the GitHub repository. Modify the input folder path to your input folder.
        `python ProteinCartography/query_uniprot.py -i uniprot_ids.txt -o input_ff/uniprot_features.tsv`
    **3.2 Manually generating the file.**
    - For each protid in the dataset, you should gather relevant protein metadata.
    - The first column of the file should be the unique protein identifer, with `protid` as the column name.
    - You should have additional columns for relevant metadata as described in [Feature file main columns](#feature-file-main-columns).
    - Default columns that are missing will be ignored.
4. Run this command, replacing the `config_ff.yml` with your config file and `8` with the number of cores you want to allocate to Snakemake.
```
snakemake --snakefile Snakefile_ff --configfile config_ff.yml --use-conda --cores 8
```

---
## Directory Structure
- [Snakefile](Snakefile): the Snakemake pipeline that orchestrates this repo's functions.
- [config.yml](config.yml): default config file for the pipeline.
- [envs/](envs/): the conda environments used for this repo.
- [demo/](demo/): contains `config_actin.yml` as well as a FASTA and PDB for human actin, used for the demo.
- [ProteinCartography/](ProteinCartography/): scripts that are called by Snakemake, importable in Python.
- [pub/](pub/): files related to the ProteinCartography pub.

---
## Pipeline Overview
This repo cotains a Snakemake pipeline that takes input `.fasta` and `.pdb` files of interest.
The rulegraph for this pipeline is as follows:  
![rulegraph](rulegraph.png)

The steps of the pipeline have the following functionality:

### Protein Folding
0. Fold input FASTA files using ESMFold.
    - The pipeline starts with input FASTA and/or PDB files, one entry per file, with a unique protein identifier (`protid`).
    - If a matching PDB file is not provided, the pipeline will use the ESMFold API to generate a PDB, if the FASTA is less than 400aa in length.

### Protein Search
1. Search the Alphafold databases using queries to the FoldSeek webserver API for each provided `.pdb` file.  
    - Currently, we're limited to a maximum of 1000 hits per protein.  

2. Search the non-redundant GenBank/RefSeq database using blastp for each provided `.fasta` file.  
    - Takes the resulting output hits and maps each GenBank/RefSeq hit to a Uniprot ID using `requests` and [the Uniprot REST API](https://rest.uniprot.org/docs/?urls.primaryName=idmapping#/job/submitJob). 
    - This can sometimes fail for unknown reasons.
    
### Download Data
3. Aggregate the list of Foldseek and BLAST hits from all input files into a single list of Uniprot IDs.  

4. Download a `.pdb` file from each protein from AlphaFold.  
    - This part is a bit slow, as it's currently limited to the number of cores provided to Snakemake.

5. Download annotation and feature information for each protein from Uniprot.  
    
### Clustering

6. Generate a similarity matrix and cluster all protein .pdb files using FoldSeek.  

7. Perform dimensionality reduction and clustering on the similarity matrix.
    - By default, we perform 30-component PCA and pass this to both TSNE and UMAP for visualization.
    - For clustering, we use the defaults of [`scanpy`'s Leiden clustering implementation](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.leiden.html#scanpy-tl-leiden).
    
### Data Analysis and Aggregation

8. Generate a variety of `_features.tsv` files.
    - Each file has, as its first column, a list of protein ids (protid) that are shared between files.
    - We query Uniprot to get all metadata from that service as a `uniprot_features.tsv` file.
    - Foldseek generates a `struclusters_features.tsv` file.
    - We perform Leiden clustering to generate a `leiden_features.tsv` file.
    - We extract from Foldseek's all-v-all TMscore analysis a distance from every protid to our input protids as `<input_protid>_distance_features.tsv` files.
    - We extract from Foldseek search a fraction sequence identy for every protid in our input protids as `<input_protid>_fident_features.tsv` files.
    - We subtract the fraction sequence identity from the TMscore to generate a `<input_protid>_convergence_features.tsv` file.
    - We determine the source of each file in the analysis (whether it was found from blast or foldseek) as the `source_features.tsv` file.
    
9. Aggregate features.
    - All of the features.tsv files are combined into one large `aggregated_features.tsv` file.

### Plotting

10. Calculate per-cluster structural similarities.
    - For every cluster of proteins, get the mean sequence similarity within that cluster and between that cluster and every other cluster.
    - Plot it as a heatmap with suffix `_leiden_similarity.html`.

11. Perform simple semantic analysis on Uniprot annotations.
    - For the annotations in each cluster, aggregate them and count the frequency of every full annotation string.
    - Perform a word count analysis on all of the annotations and generate a word cloud.
    - Save as a PDF file with suffix `_semantic_analysis.pdf`.
    
12. Build an explorable HTML visualization using `Plotly` based on the aggregated features.  
    - An example can be found [here](examples/scatter.html)
    - Each point has hover-over information
    - Default parameters include:
        - **Leiden Cluster:** Protein cluster as called by [scanpy's implementation of Leiden clustering](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html).
        - **Annotation Score:** [Uniprot annotation score](https://www.uniprot.org/help/annotation_score) from 0 to 5 (5 being the best evidence).
        - **Broad Taxon:** the broad taxonomic category for the source organism. There are two modes: 'euk' and 'bac'.
            - Assigns the "smallest" taxonomic scope from the rankings below for each point. So, a mouse gene would get `Mammalia` but not `Vertebrata`.
            - For 'euk', uses the taxonomic groups `Mammalia, Vertebrata, Arthropoda, Ecdysozoa, Lophotrochozoa, Metazoa, Fungi, Viridiplantae, Sar, Excavata, Amoebazoa, Eukaryota, Bacteria, Archaea, Viruses`
            - For 'bac', uses the taxonomic groups `Pseudomonadota, Nitrospirae, Acidobacteria, Bacillota, Spirochaetes, Cyanobacteria, Actinomycetota, Deinococcota, Bacteria, Archaea, Viruses, Metazoa, Fungi, Viridiplantae, Eukaryota`
        - **Length:** length of the protein in amino acids.
        - **Source:** how the protein was added to the clustering space (blast, foldseek or both).
        
    - Power users can customize the plots using a variety of rules, described below.

---
## Plotting Rules for `plot_interactive()`

The `plot_interactive()` function has two required arguments:
- a path to the aggregated features file with coordinates that you'd like to plot
- a `plotting_rules` dictionary describing how the data should be plotted

The `plotting_rules` dictionary should have the following format.
Each column is an entry in the dictionary containing a dictionary of rules.
```
{
    'column1.name': {
        'type': 'categorical',
        'parameter1': value,
        'parameter2': value,
        ...
    }
    'column2.name': {
        'type': 'hovertext',
        ...
    }
}
```
The possible rules for each column are as follows:
### For any plot type
- **'type'(required):**
    - `'categorical'`, `'continuous'`, `'taxonomic'`, or `'hovertext'`
    - Determines the plotting style of the data.
    - If 'categorical', expects the data to be in string format.
    - If 'continuous', expects the data to be in numerical (int or float) format.
    - If 'taxonomic', expects an ordered list of taxa to be used for grouping.
    - If 'hovertext', won't plot the data as a dropdown menu option, but will include the data in the hover-over text box.
- **'fillna':**
    - A value to fill rows that are `np.nan`.
    - For categorical plots, usually empty string `''`.
    - For continuous plots, usually `0`.
- **'apply':**
    - A function that will be applied to every element of the feature before plotting.
    - This can be used to convert a numerical value into a categorical value with a lambda function.
    - For example, `lambda x: str(x)`
- **'skip_hover':**
    - Boolean, whether to include this data in the hovertext.
- **'textlabel':**
    - A string that replaces the column name on buttons and in hover text.
    - Useful if the column name is too ugly.

### For 'continuous' plots
- **'color_scale':**
    - A [named Plotly colorscale](https://plotly.com/python/builtin-colorscales/), or:
    - A [list of lists that can be converted into a Plotly colorscale](https://plotly.com/python/colorscales/#explicitly-constructing-a-color-scale) (don't use tuples).
- **'cmin':**
    - Minimum value to use for color scale. Defaults to minimum value for that column.
- **'cmax':**
    - Maximum value to use for color scale. Defaults to maximum value for that column.
- Note: if the value of 'fillna' is lower than the cmin, na values will have their own light grey coloration, distinct from the rest of the color scale. 
  Good for indicating values that are missing, rather than actually calculated.

### For 'categorical' and 'taxonomic' plots
- **'color_order':**
    - A list of HEX code colors used for coloring categorical variables.
    - If there aren't enough colors for unique values of the data, will generate up to 3x more colors of varying darkness.
- **'color_dict':**
    - A dictionary of key: value pairs for coloring categorical variables.
    - The key is the name of the category and the value is a HEX code color.

### For 'taxonomic' plots
- **'taxon_order':**
    - Exclusively used for 'taxonomic' style plots. A list of ranked-order taxa for categorization.
    - Should start with more-specific taxa and expand to less-specific taxa.

---
## File Conventions

The pipeline generates a large number of .txt and .tsv files with specific formatting expectations.  
Many of the pipeline's scripts accept these specific format conventions as input or return them as output.  
These are the primary formats and their descriptions.

### Accession list files (ACC)
These files end with `'.txt'` and contain a list of accessions (RefSeq, GenBank, Uniprot), one per line.

- **Example:**
    ```
    A0A2J8L4A7
    K7EV54
    A0A2J8WJR8
    A0A811ZNA7
    ...
    ```
- **Input to:** 
    - [`aggregate_lists.py`](ProteinCartography/aggregate_lists.py)
    - [`get_source.py`](ProteinCartography/get_source.py)
    - [`map_refseqids.py`](ProteinCartography/map_refseqids.py)
    - [`query_uniprot.py`](ProteinCartography/query_uniprot.py)
    - [`rescue_mapping.py`](ProteinCartography/rescue_mapping.py)
- **Output from:**
    - [`aggregate_lists.py`](ProteinCartography/aggregate_lists.py)
    - [`extract_blasthits.py`](ProteinCartography/extract_blasthits.py)
    - [`extract_foldseekhits.py`](ProteinCartography/extract_foldseekhits.py)
    - [`map_refseqids.py`](ProteinCartography/map_refseqids.py)
    - [`rescue_mapping.py`](ProteinCartography/rescue_mapping.py)

### Matrix File (MTX)
These files end with `'.tsv'` and contain distance or similarity matrices, usually all-v-all.  

- **Example:**:
    | protid | A0A2J8L4A7 | K7EV54 | A0A2J8WJR8 | A0A811ZNA7 | 
    |-------:|:----------:|:------:|:----------:|:----------:|
    | A0A2J8L4A7 |  1   | 0.9  | 0.85 | 0.7  |
    | K7EV54     | 0.9  |  1   | 0.91 | 0.6  |
    | A0A2J8WJR8 | 0.85 | 0.91 |  1   | 0.71 |
    | A0A811ZNA7 | 0.7  | 0.6  | 0.71 |  1   |
- **Input to:** 
    - [`cluster_similarity.py`](ProteinCartography/cluster_similarity.py)
    - [`dim_reduction.py`](ProteinCartography/dim_reduction.py)
    - [`extract_input_distances.py`](ProteinCartography/extract_input_distances.py)
- **Output from:**
    - [`foldseek_clustering.py`](ProteinCartography/foldseek_clustering.py)

### Features files (FTF)
These files end with `'.tsv'` and contain a `protid` column, which is the unique identifier of each protein in the dataset.  
The remaining columns are annotations for each protein. These annotations can be of any data type.

- **Example:**
    | protid | Length | LeidenCluster | Organism |
    |-------:|:---------------:|:-------------:|:-----------------------:|
    | A0A2J8L4A7 | 707 | 1 | Pan troglodytes (Chimpanzee) |
    | K7EV54 | 784 | 2 | Pongo abelii (Sumatran orangutan) |
    | A0A2J8WJR8 | 707 | 1 | Pongo abelii (Sumatran orangutan) |
    | A0A811ZNA7 | 781 | 1 | Nyctereutes procyonoides (Raccoon dog) |
- **Input to:**
    - [`aggregate_features.py`](ProteinCartography/aggregate_features.py)
    - [`calculate_convergence.py`](ProteinCartography/calculate_convergence.py)
    - [`dim_reduction.py`](ProteinCartography/dim_reduction.py) 
    - [`plot_interactive.py`](ProteinCartography/plot_interactive.py)
    - [`semantic_analysis.py`](ProteinCartography/semantic_analysis.py)
- **Output from:**
    - [`aggregate_fident.py`](ProteinCartography/aggregate_fident.py)
    - [`calculate_convergence.py`](ProteinCartography/calculate_convergence.py)
    - [`extract_input_distances.py`](ProteinCartography/extract_input_distances.py)
    - [`foldseek_clustering.py`](ProteinCartography/foldseek_clustering.py)
    - [`get_source.py`](ProteinCartography/get_source.py)
    - [`leiden_clustering.py`](ProteinCartography/leiden_clustering.py)
    - [`query_uniprot.py`](ProteinCartography/query_uniprot.py)

#### Feature file main columns

A variety of metadata features for each protein are usually pulled from Uniprot for visualization purposes.  
The following table describes those features with examples.  

If you are providing a set of custom proteins (such as those not fetched from Uniprot), you may want to include a `features_override.tsv` file that contains these features for your proteins of interest. This will allow you to visualize your protein correctly in the interactive HTML map.  

Features used for color schemes in the default plotting rules are marked with *(Plotting)* below. Features used only for hover-over description are marked with *(Hovertext)*.

| feature | example | description | source |
|--------:|:-------:|:------------|:-------|
|`"protid"` | `"P42212"` | *(Required)* the unique identifier of the protein. Usually the Uniprot accession, but can be any alphanumeric string | User-provided or Uniprot |
|`"Protein names"`| `"Green fluorescent protein"` | *(Hovertext)* a human-readable description of the protein | Uniprot |
|`"Gene Names (primary)"` | `"GFP"` | *(Hovertext)* a gene symbol for the protein | Uniprot |
|`"Annotation"`| `5` | *(Plotting)* Uniprot [Annotation Score](https://www.uniprot.org/help/annotation_score) (0 to 5) | Uniprot |
|`"Organism"` | `"Aequorea victoria (Water jellyfish) (Mesonema victoria)"` | *(Hovertext)* Scientific name (common name) (synonyms) | Uniprot |
|`"Taxonomic lineage"`|`"cellular organisms (no rank), Eukaryota (superkingdom), ... Aequoreidae (family), Aequorea (genus)"`| string of comma-separated `Lineage name (rank)` for the organism's full taxonomic lineage | Uniprot |
|`"Lineage"`|`["cellular organisms", "Eukaryota", ... "Aequoreidae", "Aequorea"]` | *(Plotting)* ordered list of lineage identifiers without rank information, generated from `"Taxonomic lineage"` | ProteinCartography |
|`"Length"`| `238` | *(Plotting)* number of amino acids in protein | Uniprot |
