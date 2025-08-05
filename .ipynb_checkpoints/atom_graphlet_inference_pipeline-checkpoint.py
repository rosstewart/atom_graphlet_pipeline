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
from postprocessing_utils import get_significant_predictions
import subprocess
import os
import sys
import argparse

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

property_help = "Property codes:\n" + "\n".join([f"  {k}: {v}" for k, v in property_descs.items()])

parser = argparse.ArgumentParser(
    description='Predict residue properties',
    epilog=property_help,
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument('--input', required=True, type=str,
                    help='Input PDB file')
parser.add_argument('--output', required=True, type=str,
                    help='Output numpy file for predictions')
parser.add_argument('--property', required=True, type=str,
                    choices=property_descs.keys(),
                    help='Property code (see list below)')

parser.add_argument('--p_threshold', type=float, default=0.05,
                    help='P-value threshold for significance (default: 0.05)')
parser.add_argument('--chain', type=str, default='A',
                    help='PDB chain to analyze (default: A)')

args = parser.parse_args()

if not os.path.exists(args.input):
    print(f"Error: Input file '{args.input}' not found")
    sys.exit(1)

print(f"Processing {args.input} for {property_descs[args.property]} (chain {args.chain})...")


pdb_f = args.input
preds_f = args.output
property_code = args.property
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
    # backbone atoms
    "N": "N", "CA": "A", "C": "C", "O": "O",
    # beta carbon
    "CB": "B",
    # gamma atoms
    "CG": "G", "OG": "G", "SG": "G",
    # delta atoms  
    "CD": "D", "OD": "D", "SD": "D", "ND": "D",
    # epsilon atoms
    "CE": "E", "OE": "E", "NE": "E",
    # zeta atoms
    "CZ": "Z", "NZ": "Z",
    # eta atoms
    "NH": "H",
    # special cases
    "OH": "H",  # tyrosine hydroxyl
}



parser = PDB.PDBParser(QUIET=True)
wd = os.getcwd()
graphlet_wd = './graphlet_counting'
model_dir = './svm'
model_path = f'{model_dir}/{property_code}_calibrated_svm.pkl'
cv_results_dir='./crossval_results'
edge_dist_threshold = 7.5

with tempfile.TemporaryDirectory() as save_dir:
    
    # generate contact graph
    pdb_id, mat_data = make_graph(pdb_f, edge_dist_threshold, parser, three_letter_to_one, pdb_atom_mapping, chain=args.chain)
    
    # write files needed downstream
    write_graph_and_labels(pdb_id, mat_data, save_dir)
    make_pos_file(pdb_id, save_dir)

    # generate graphlet features (7.5 Å edge dist threshold, only up to 4-graphlet for speed)
    os.chdir(graphlet_wd) # need to be in graphlet dir unless i want to refactor more stuff
    subprocess.run([f'./run_atom_std.sh', pdb_id, save_dir], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir(wd)
    
    svml_to_sparse(save_dir, graphlet_idx_mapping_f=f'{model_dir}/{property_code}_graphlet_idx_mapping.pkl')
    residue_preds = run_svm_inference(pdb_id, save_dir, model_path)

    empirical_pvalues, significant_indices = get_significant_predictions(property_code, residue_preds, cv_results_dir=cv_results_dir, p_threshold=args.p_threshold)

    np.save(preds_f, residue_preds)

print(f'\nresidue predictions saved to {preds_f}\n')
