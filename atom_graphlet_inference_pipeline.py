'''
Ross Stewart 05/30/25

ncac_cat_pipeline

predicts residue property probabilities from an input pdb file

args:
figure it out

output:
posteriors
'''


import numpy as np
from Bio import PDB
import tempfile
from contact_graph_utils import make_graph, write_graph_and_labels, make_pos_file
from svm_utils import svml_to_sparse, run_svm_inference
import subprocess
import os
import sys

property_descs = {
    'ADP': 'ADP-binding',
    'Allo': 'Allosteric residue',
    'ATP': 'ATP-binding',
    'CA': 'Calcium-binding',
    'Cat': 'Catalytic residue',
    'CD': 'Cadmium-binding',
    'CO': 'Cobalt-binding',
    'CU': 'Copper-binding',
    'DNA': 'DNA-binding',
    'FAD': 'FAD-binding',
    'FE': 'Iron-binding',
    'FMN': 'FMN-binding',
    'GDP': 'GDP-binding',
    'GTP': 'GTP-binding',
    'HEM': 'Heme-binding',
    'Hotspot': 'Protein-protein interaction hotspot',
    'K': 'Potassium-binding',
    'MG': 'Magnesium-binding',
    'MN': 'Manganese-binding',
    'NA': 'Sodium-binding',
    'NAD': 'NAD-binding',
    'Nglyco': 'N-linked glycosylation site',
    'NI': 'Nickel-binding',
    'Phos': 'Phosphorylation site',
    'PLP': 'PLP-binding',
    'PPI': 'Protein-protein interaction site',
    'RNA': 'RNA-binding',
    'UDP': 'UDP-binding',
    'ZN': 'Zinc-binding',
 }

if len(sys.argv) != 4:
    print(f'\nUsage: python {sys.argv[0]} <input.pdb> <output.npy> <property_code>\n\nPlease choose one of the property codes below:\n{property_descs}\n')
    sys.exit(1)



pdb_f = sys.argv[1]
assert os.path.exists(pdb_f)
preds_f = sys.argv[2]
property_code = sys.argv[3]
assert property_code in property_descs

three_letter_to_one = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Asx": "B", "Glx": "Z",
    "Xaa": "X", "Ter": "*"
}

pdb_atom_mapping = {
    # Backbone atoms
    "N": "N", "CA": "A", "C": "C", "O": "O",
    # Beta carbon
    "CB": "B",
    # Gamma atoms
    "CG": "G", "OG": "G", "SG": "G",
    # Delta atoms  
    "CD": "D", "OD": "D", "SD": "D", "ND": "D",
    # Epsilon atoms
    "CE": "E", "OE": "E", "NE": "E",
    # Zeta atoms
    "CZ": "Z", "NZ": "Z",
    # Eta atoms
    "NH": "H",
    # Special cases
    "OH": "H",  # Tyrosine hydroxyl
}



parser = PDB.PDBParser(QUIET=True)
wd = os.getcwd()
graphlet_wd = './graphlet_counting'
model_path = f'./svm/{property_code}_calibrated_svm.pkl'
edge_dist_threshold = 7.5

with tempfile.TemporaryDirectory() as save_dir:
    
    # generate contact graph
    pdb_id, mat_data = make_graph(pdb_f, edge_dist_threshold, parser, three_letter_to_one, pdb_atom_mapping)
    
    # write files needed downstream
    write_graph_and_labels(pdb_id, mat_data, save_dir)
    make_pos_file(pdb_id, save_dir)

    # generate graphlet features (7.5 Å edge dist threshold, only up to 4-graphlet for speed)
    os.chdir(graphlet_wd) # need to be in graphlet dir unless i want to refactor more stuff
    subprocess.run([f'./run_atom_std.sh', pdb_id, save_dir], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir(wd)
    
    svml_to_sparse(save_dir, graphlet_idx_mapping_f=f'{model_dir}/{property_code}_graphlet_idx_mapping.pkl')
    residue_preds = run_svm_inference(pdb_id, save_dir, model_path)

    np.save(preds_f, residue_preds)

print(f'\n residue predictions saved to {preds_f}\n')
