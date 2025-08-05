# Atom Graphlet Property Predictor Pipeline
Pipeline to predict up to 29 residue properties of a protein using atom graphlets

Usage: `python atom_graphlet_inference_pipeline.py <input.pdb> <output.npy> <property_code>`

Choose from one of the property codes below (each mapped to the text description and cross-validation PU AUC)

    'ADP': 'ADP-binding (AUC: 0.860)',
    'Allo': 'Allosteric residue (AUC: 0.693)',
    'ATP': 'ATP-binding (AUC: 0.834)',
    'CA': 'Calcium-binding (AUC: 0.876)',
    'Cat': 'Catalytic residue (AUC: 0.893)',
    'CD': 'Cadmium-binding (AUC: 0.850)',
    'CO': 'Cobalt-binding (AUC: 0.899)',
    'CU': 'Copper-binding (AUC: 0.937)',
    'DNA': 'DNA-binding (AUC: 0.828)',
    'FAD': 'FAD-binding (AUC: 0.883)',
    'FE': 'Iron-binding (AUC: 0.952)',
    'FMN': 'FMN-binding (AUC: 0.820)',
    'GDP': 'GDP-binding (AUC: 0.869)',
    'GTP': 'GTP-binding (AUC: 0.673)',
    'HEM': 'Heme-binding (AUC: 0.901)',
    'Hotspot': 'Protein-protein interaction hotspot (AUC: 0.756)',
    'K': 'Potassium-binding (AUC: 0.811)',
    'MG': 'Magnesium-binding (AUC: 0.836)',
    'MN': 'Manganese-binding (AUC: 0.915)',
    'NA': 'Sodium-binding (AUC: 0.747)',
    'NAD': 'NAD-binding (AUC: 0.890)',
    'Nglyco': 'N-linked glycosylation site (AUC: 0.821)',
    'NI': 'Nickel-binding (AUC: 0.917)',
    'Phos': 'Phosphorylation site (AUC: 0.709)',
    'PLP': 'PLP-binding (AUC: 0.940)',
    'PPI': 'Protein-protein interaction site (AUC: 0.792)',
    'RNA': 'RNA-binding (AUC: 0.793)',
    'UDP': 'UDP-binding (AUC: 0.875)',
    'ZN': 'Zinc-binding (AUC: 0.921)',


