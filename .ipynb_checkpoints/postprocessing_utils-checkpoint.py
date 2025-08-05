import numpy as np

def convert_to_posterior_pu(pred, alpha, p_s0, p_s1):
    """
    Corrected PU calibration using Elkan & Noto approach.
    """
    # calculate labeling frequency
    c = p_s1 / (p_s1 + alpha * p_s0)
    
    # simple calibration
    posterior = pred / c
    
    return np.clip(posterior, 0, 1)

def get_significant_predictions(prefix, new_predictions, cv_results_dir, p_threshold=0.05):
    """
    Calculate empirical p-values for new predictions based on the unlabeled distribution
    from cross-validation results, and identify significant predictions.
    
    Parameters:
    - prefix: the dataset prefix (e.g., 'Cat', 'ATP', etc.)
    - new_predictions: numpy array of new prediction scores
    - p_threshold: p-value threshold (default 0.05)
    - cv_results_dir: directory containing crossval results
    
    Returns:
    - empirical_pvalues: array of p-values for each prediction
    - significant_indices: 1-based indices of significant predictions
    """

    class_prior_dict = {
       'CA': 0.83,
       'CD': 1.32,
       'CO': 0.58,
       'CU': 0.40,
       'FE': 0.17,
       'K': 1.20,
       'MG': 0.40,
       'MN': 0.01,
       'NA': 3.62,
       'NI': 0.83,
       'ZN': 1.32,
       'Nglyco': 13.3,
       'Phos': 0.44,
       'Cat': 3.16,
       'DNA': 0.48,
       'RNA': 0.36,
       'PPI': 16.3,
       'Hotspot': 0.50,
       'ADP': 2.09,
       'ATP': 2.95,
       'FAD': 2.31,
       'FMN': 0.45,
       'GDP': 0.83,
       'GTP': 1.22,
       'HEM': 2.76,
       'NAD': 2.06,
       'PLP': 0.72,
       'UDP': 4.25,
       'Allo': 4.85,
    }
    
    # load the crossval results to get the null distribution
    labels_file = f"{cv_results_dir}/{prefix}_7_5_atom_avg_4graphlet_svm_labels.npy"
    preds_file = f"{cv_results_dir}/{prefix}_7_5_atom_avg_4graphlet_svm_preds.npy"
    
    try:
        labels = np.load(labels_file)
        preds = np.load(preds_file)
    except FileNotFoundError:
        print(f"Error: Could not find files for prefix '{prefix}'")
        return None, None
    
    # get unlabeled predictions as null distribution
    unlabeled_preds = preds[labels == -1]
    s1_size, s0_size = len(preds[labels == 1]), len(unlabeled_preds)
    
    if len(unlabeled_preds) == 0:
        print(f"Error: No unlabeled samples found for prefix '{prefix}'")
        return None, None

    # calibrate for PU setting
    unlabeled_preds = np.array([convert_to_posterior_pu(pred, class_prior_dict[prefix]/100, s0_size, s1_size) for pred in unlabeled_preds])
    new_predictions = np.array([convert_to_posterior_pu(pred, class_prior_dict[prefix]/100, s0_size, s1_size) for pred in new_predictions])
    
    # calculate empirical p-values
    # p-value = fraction of unlabeled predictions >= each new prediction
    empirical_pvalues = np.zeros(len(new_predictions))
    for i, pred in enumerate(new_predictions):
        empirical_pvalues[i] = np.mean(unlabeled_preds >= pred)
    
    # find significant predictions (p-value < threshold)
    significant_mask = empirical_pvalues < p_threshold
    significant_indices = np.where(significant_mask)[0] + 1  # Convert to 1-based
    
    # print results
    print(f"significant predictions at p-value threshold: {p_threshold}: {len(significant_indices)} out of {len(new_predictions)}")
    
    if len(significant_indices) > 0:
        print(f"\nsignificant prediction indices (1-based):")
        for idx in significant_indices:
            print(f"  residue {idx}: score={new_predictions[idx-1]:.4f}, p-value={empirical_pvalues[idx-1]:.4f}")
    else:
        print("\nno predictions meet the significance threshold")
    
    # also print the score threshold corresponding to the p-value
    score_threshold = np.percentile(unlabeled_preds, (1 - p_threshold) * 100)
    print(f"\nscore threshold for p={p_threshold}: {score_threshold:.4f}")
    print(f"(predictions above {score_threshold:.4f} are significant at p<{p_threshold})")
    
    return new_predictions, empirical_pvalues, significant_indices
