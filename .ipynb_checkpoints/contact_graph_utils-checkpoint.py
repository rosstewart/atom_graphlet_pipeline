import numpy as np
import scipy.sparse as sp
from scipy.io import savemat # For matlab!
import os
from Bio import PDB
import glob
import tempfile
import re
from scipy.io import loadmat

def map_pdb_atom(atom_name, pdb_atom_mapping):
    """Maps a PDB atom code to a single letter representation, ignoring hydrogens"""
    # Skip hydrogens entirely
    if atom_name.startswith('H') or atom_name[0].isdigit() and 'H' in atom_name:
        return None
    
    # Remove numbers (e.g., CG1 → CG, OD1 → OD)
    base_atom = re.sub(r'\d+', '', atom_name)
    
    # Return mapped value or "?" if not found
    return pdb_atom_mapping.get(base_atom, "?")

'''
make atom contact graph and extract necessary information
'''
def make_graph(pdb_f, edge_dist_threshold, parser, three_letter_to_one, pdb_atom_mapping):
    '''
    construct graph with one chain in one graph
    '''

    pdb_id = pdb_f.split('/')[-1].split('.')[0]

    model = parser.get_structure(pdb_id, pdb_f)[0]['A']
    
    '''
    get residues and coords
    '''
    # aa_labels = []
    atom_labels = []
    atom_coords = []
    atom_res_indices = []
    num_residues = 0
    num_atoms = 0
    light_atoms = ('H')
    for residue in model:
        if PDB.is_aa(residue):
            num_residues += 1
            # if residue.get_resname().capitalize() in three_letter_to_one:
            #     aa_labels.append(three_letter_to_one[residue.get_resname().capitalize()])
            
            for atom in residue:
                if map_pdb_atom(atom.get_name(), pdb_atom_mapping) is not None: # is not hydrogen
                    atom_labels.append(map_pdb_atom(atom.get_name(), pdb_atom_mapping))
                    atom_coords.append(residue[atom.get_name()].coord)
                    atom_res_indices.append(num_residues - 1)
                    num_atoms += 1
    
    '''
    construct atom edge matrix, no self-loops for graphlet input
    '''
    edge_mat = np.zeros((num_atoms,num_atoms)) # dont include diagonal for pedja's profile kernel input
    for i in range(num_atoms):
        if atom_coords[i] is None:
            print(pdb_id,'residue',atom_res_indices[i],'atom',atom_labels[i],'has no coords')
            raise ValueError
        for j in range(i+1,num_atoms):
            if atom_coords[j] is None:
                print(pdb_id,'residue',atom_res_indices[j],'atom',atom_labels[j],'has no coords')
                raise ValueError
            
            dist = np.linalg.norm(atom_coords[i] - atom_coords[j])
            if dist <= edge_dist_threshold: 
                assert edge_mat[i,j] == edge_mat[j,i] and edge_mat[i,j] == 0
                edge_mat[i,j],edge_mat[j,i] = 1,1 # undirected edges
    
    '''
    save to matlab input file
    assume matlab automatically converts zero-based indexing to one-based
    '''
    # savemat(f'{save_dir}/{pdb_id}.mat', {'G': sp.csr_matrix(edge_mat), 'L': aa_labels, 'AL': atom_labels, 'ARI': atom_res_indices}) # aa labels is the same as pdb seq as a list
    # mat = {'G': sp.csr_matrix(edge_mat), 'L': aa_labels, 'AL': atom_labels, 'ARI': atom_res_indices}
    mat = {'G': sp.csr_matrix(edge_mat), 'AL': atom_labels, 'ARI': atom_res_indices}
    return pdb_id, mat

'''
write input files for graphlet counting script
'''
def write_graph_and_labels(pdb_id, mat_data, save_dir):
    dense_mat = mat_data['G'].toarray()
    # pdb_seq = mat_data['L']
    atom_seq = mat_data['AL']
    atom_res_indices = mat_data['ARI']
    
    # if len(pdb_seq) == 1:
    #     pdb_seq = pdb_seq[0]
    # else:
    #     pdb_seq = ''.join(pdb_seq)
    
    if len(atom_seq) == 1:
        atom_seq = atom_seq[0]
    else:
        atom_seq = ''.join(atom_seq)
    
    if len(atom_res_indices) == 1:
        atom_res_indices = atom_res_indices[0]
    
    with open(f"{save_dir}/{pdb_id}.graph",'w') as f_w:
        for i in range(len(dense_mat)):
            nonzero = dense_mat[i].nonzero()
            assert len(nonzero) == 1
            close_indices = '\t'.join(map(str, nonzero[0]))
            f_w.write(f"{i} {close_indices}\n")
    
    # with open(f"{save_dir}/{pdb_id}.aa_labels",'w') as f_w:
    #     f_w.write(pdb_seq)
    with open(f"{save_dir}/{pdb_id}.atom_labels",'w') as f_w:
        f_w.write(atom_seq)
    with open(f"{save_dir}/{pdb_id}.atom_res_indices",'w') as f_w:
        f_w.write(','.join(np.array(atom_res_indices).astype(str).tolist()))

'''
mark all atoms as positive for inference to ensure they are predicted on
'''
def make_pos_file(pdb_id, save_dir):
    with open(f"{save_dir}/{pdb_id}.atom_pos",'w') as f_w:
        with open(f"{save_dir}/{pdb_id}.atom_labels",'r') as f:
            for atom_idx in range(len(f.readlines()[0])):
                f_w.write(f'{atom_idx}\n')

