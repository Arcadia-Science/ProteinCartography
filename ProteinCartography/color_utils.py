import arcadia_pycolor as apc

##############################
## Standard plotting colors ##
##############################

arcadia_viridis = apc.Gradients["arcadia:viridis"].grad_nested_list
arcadia_viridis_r = apc.Gradients["arcadia:viridis_r"].grad_nested_list

arcadia_magma = apc.Gradients["arcadia:magma"].grad_nested_list
arcadia_magma_r = apc.Gradients["arcadia:magma_r"].grad_nested_list

arcadia_cividis = apc.Gradients["arcadia:cividis"].grad_nested_list
arcadia_cividis_r = apc.Gradients["arcadia:cividis_r"].grad_nested_list

arcadia_poppies = apc.Gradients["arcadia:poppies"].grad_nested_list
arcadia_poppies_r = apc.Gradients["arcadia:poppies_r"].grad_nested_list

arcadia_pansies = apc.Gradients["arcadia:pansies"].grad_nested_list
arcadia_pansies_r = apc.Gradients["arcadia:pansies_r"].grad_nested_list

arcadia_dahlias = apc.Gradients["arcadia:dahlias"].grad_nested_list
arcadia_dahlias_r = apc.Gradients["arcadia:dahlias_r"].grad_nested_list

plddt_gradient_dict = {
    "color_dict": apc.dragon
    | apc.amber
    | apc.canary
    | apc.vitalblue
    | {"arcadia:cobalt": "#4A72B0"},
    "values": [0, 0.25, 0.6, 0.8, 1],
}

# instantiate a new Gradient object
PLDDT_GRADIENT = apc.Gradient(
    name="my_gradient",
    color_dict=plddt_gradient_dict["color_dict"],
    values=plddt_gradient_dict["values"],
)

PLDDT_CMAP = PLDDT_GRADIENT.grad_nested_list

RESIDUE_CONFIDENCE_COLORS = {
    "very_high": "#4A72B0",
    "confident": apc.All["arcadia:vitalblue"],
    "low": apc.All["arcadia:canary"],
    "very_low": apc.All["arcadia:amber"],
}

ANNOTATION_SCORE_COLORS = [
    apc.All["arcadia:brightgrey"],
    apc.All["arcadia:aster"],
    apc.All["arcadia:aegean"],
    apc.All["arcadia:seaweed"],
    apc.All["arcadia:lime"],
    apc.All["arcadia:canary"],
]

ANNOTATION_SCORE_COLOR_DICT = dict(zip([str(i) for i in range(6)], ANNOTATION_SCORE_COLORS))

EUK_COLOR_DICT = {
    "Mammalia": apc.All["arcadia:oat"],
    "Vertebrata": apc.All["arcadia:canary"],
    "Arthropoda": apc.All["arcadia:seaweed"],
    "Ecdysozoa": apc.All["arcadia:mint"],
    "Lophotrochozoa": apc.All["arcadia:aegean"],
    "Metazoa": apc.All["arcadia:amber"],
    "Fungi": apc.All["arcadia:chateau"],
    "Viridiplantae": apc.All["arcadia:lime"],
    "Sar": apc.All["arcadia:rose"],
    "Excavata": apc.All["arcadia:wish"],
    "Amoebazoa": apc.All["arcadia:periwinkle"],
    "Eukaryota": apc.All["arcadia:aster"],
    "Bacteria": apc.All["arcadia:slate"],
    "Archaea": apc.All["arcadia:dragon"],
    "Viruses": apc.All["arcadia:denim"],
}

BAC_COLOR_DICT = {
    "Pseudomonadota": apc.All["arcadia:periwinkle"],
    "Nitrospirae": apc.All["arcadia:vitalblue"],
    "Acidobacteria": apc.All["arcadia:mars"],
    "Bacillota": apc.All["arcadia:mint"],
    "Spirochaetes": apc.All["arcadia:aegean"],
    "Cyanobacteria": apc.All["arcadia:seaweed"],
    "Actinomycetota": apc.All["arcadia:canary"],
    "Deinococcota": apc.All["arcadia:rose"],
    "Bacteria": apc.All["arcadia:slate"],
    "Archaea": apc.All["arcadia:dragon"],
    "Viruses": apc.All["arcadia:denim"],
    "Metazoa": apc.All["arcadia:amber"],
    "Fungi": apc.All["arcadia:chateau"],
    "Viridiplantae": apc.All["arcadia:lime"],
    "Eukaryota": apc.All["arcadia:wish"],
}

SOURCE_COLOR_DICT = {
    "blast": apc.All["arcadia:canary"],
    "foldseek": apc.All["arcadia:aegean"],
    "blast+foldseek": apc.All["arcadia:amber"],
    "None": apc.All["arcadia:brightgrey"],
}

PDB_ORIGIN_COLOR_DICT = {
    "AlphaFold": "#4A72B0",
    "ESMFold": apc.All["arcadia:vitalblue"],
    "PDB": apc.All["arcadia:bluesky"],
    "Other": apc.All["arcadia:marineblue"],
    "None": apc.All["arcadia:brightgrey"],
}
