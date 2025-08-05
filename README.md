# Atom Graphlet Property Predictor Pipeline
Pipeline to predict up to 29 residue properties of a protein using atom graphlet SVM and give empirical p-values for predictions based on the unlabeled dataset

## Usage

`python atom_graphlet_inference_pipeline.py --input <pdb_file> --output <output_file> --property <property_code> [options]`

### Required Arguments

`--input` : Path to input PDB file

The protein structure file in PDB format.
Example: `--input protein.pdb`


`--output` : Path to output numpy file

Where to save the prediction scores (numpy array).
Example: `--output predictions.npy`


`--property` : Property code to predict

Specifies which functional property to predict.
Must be one of the valid property codes (see list below).
Example: `--property Cat`

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

### Optional Arguments

`--p_threshold` : P-value threshold for significance (default: `0.05`)

Controls the statistical significance cutoff.
Lower values = more stringent filtering.
Example: `--p_threshold 0.01`


`--chain` : PDB chain to analyze (default: `'A'`)

Specifies which chain to analyze in multi-chain proteins.
Example: `--chain B`



## Examples

### Basic Usage

Predict catalytic residues in chain A with default p-value threshold:
`python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output cat_predictions.npy --property Cat`

### Custom P-value Threshold

Use a more stringent p-value cutoff of 0.01:
`python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output atp_predictions.npy --property ATP --p_threshold 0.01`

### Specific Chain

Analyze chain B instead of the default chain A:
`python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output phos_predictions.npy --property Phos --chain B`

### Multiple Properties

To predict multiple properties, run the script multiple times:

```
# Predict catalytic residues
python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output cat_predictions.npy --property Cat

# Predict zinc-binding residues
python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output zn_predictions.npy --property ZN

# Predict phosphorylation sites with stringent threshold
python atom_graphlet_inference_pipeline.py --input 1xyz.pdb --output phos_predictions.npy --property Phos --p_threshold 0.001
```

## Output

The tool generates two output files:

1. Prediction scores (`<output>.npy`):

    • Numpy array containing prediction scores for each residue

    • Higher scores indicate higher confidence in the predicted property


2. Significant residue indices (`<output>_significant_res_indices.npy`):

    • Numpy array listing residue indices that meet the p-value threshold

    • Contains 1-based residue indices

## Notes

• Residue indices in the output are 1-based (first residue = 1)

• The p-value threshold determines statistical significance based on empirical null distributions

• Lower p-values (e.g., 0.01, 0.001) provide more conservative predictions

• The tool analyzes one chain at a time; for multi-chain analysis, run separately for each chain

• Some properties (e.g. Phos, Nglyco) only occur on specific amino acids. The training sets are designed appropriately, however this program DOES NOT account for filtering based on amino acid. Please postprocess your predictions accordingly.

• You must run the script from the directory is it located, otherwise you may have to refactor to add absolute paths.

