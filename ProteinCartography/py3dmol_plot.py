#!/usr/bin/env python
import os
from Bio.PDB import CEAligner, PDBParser, Selection, cealign, PDBIO
import py3Dmol

def plot_structure_pair(file1: str, file2: str, 
                        color1 = 'red', color2 = 'blue',
                        width = 800, height = 600, style = 'cartoon',
                        label_dict = {}, keep_aligned = False):
    '''
    Interactively plots a pair of structures from pdb file paths.
    Works in Jupyter notebooks.
    
    Create a temporary reference and aligned structure in the same location as the starting proteins.
    Deletes these temporary files by default.
    
    Args:
        file1 (str): path to first pdb file (reference structure)
        file2 (str): path to second pdb file
        color1 (str): matplotlib color name or HEX code
        color2 (str): matplotlib color name or HEX code
        width (int): width of visualization in pixels
        height (int): height of visualization in pixels
        style (str): 'cartoon' or 'stick'
        label_dict (dict): key-value pairs for plot labels where key is filepath and value is label text
        keep_aligned (bool): whether to keep the aligned versions of the PDBs, default False.
    '''
    #get filepath from pdb name
    pdb1_filename = os.path.abspath(file1)
    pdb2_filename = os.path.abspath(file2)
    
    # Load the PDB files
    pdb1 = PDBParser().get_structure('pdb1', pdb1_filename)
    pdb2 = PDBParser().get_structure('pdb2', pdb2_filename)
    
    #align the structures with cealigner
    aligner=CEAligner()
    aligner.set_reference(pdb1)
    aligner.align(pdb2, transform=True)
    
    reference_temp_file1 = file1 + '.reference_structure.pdb'
    aligned_temp_file2 = file2 + '.aligned_structure.pdb'
    
    io = PDBIO()
    io.set_structure(pdb1)
    io.save(reference_temp_file1)
    io.set_structure(pdb2)
    io.save(aligned_temp_file2)

    # Initialize the view object
    view = py3Dmol.view(width=width, height=height)
    
    # Load the PDB files
    with open(reference_temp_file1, 'r') as f:
        pdb_data1 = f.read()
    with open(aligned_temp_file2, 'r') as f:
        pdb_data2 = f.read()
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data1, 'pdb')
    view.addModel(pdb_data2, 'pdb')
    
    if not keep_aligned:
        os.remove(reference_temp_file1)
        os.remove(aligned_temp_file2)
    
    if style == 'cartoon' or style == 'stick':
        view.setStyle({'model': 0}, {style: {'color': color1}})
        view.setStyle({'model': 1}, {style: {'color': color2}})
    else:
        raise Exception('Must choose "cartoon" or "stick" as style.')
        
    if label_dict != {}:
        label1 = label_dict[file1]
        label2 = label_dict[file2]
    else:
        label1 = os.path.basename(file1)
        label2 = os.path.basename(file2)
    
    view.addLabel(label1, {'backgroundColor': 'white', 
                           'backgroundOpacity': 0.1, 
                           'fontColor': color1, 
                           'useScreen': True, 
                           'position': {'x': 0, 'y': 10, 'z': 0}})
    view.addLabel(label2, {'backgroundColor': 'white', 
                           'backgroundOpacity': 0.1, 
                           'fontColor': color2, 
                           'useScreen': True, 
                           'position': {'x': 0, 'y': 30, 'z': 0}})
    view.addLabel(f"RMSD: {aligner.rms: .2f}", {'backgroundColor': 'white', 
                           'backgroundOpacity': 0.1, 
                           'fontColor': 'black', 
                           'useScreen': True, 
                           'position': {'x': 0, 'y': 60, 'z': 0}})
    view.zoomTo()
    return view.show()
