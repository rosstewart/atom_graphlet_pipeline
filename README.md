# Atom Graphlet Property Predictor Pipeline
Pipeline to predict up to 29 residue properties of a protein using atom graphlets

Usage: `python atom_graphlet_inference_pipeline.py <input.pdb> <output.npy> <property_code>`

Choose from one of the property codes below (each mapped to the text description)

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

