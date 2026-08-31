"""
Generic CSV/TSV reader with automatic column identification for virtual screening results.
"""

from typing import Tuple, Dict, Any, Optional
import os
import numpy as np
import pandas as pd


def load_docking_csv(
    filepath: str,
    score_column: Optional[str] = None,
    label_column: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, Any]]:
    """
    Loads virtual screening tabular results with automated column discovery.

    Parameters
    ----------
    filepath : str
        Path to CSV/TSV table.
    score_column : str, optional
        Name of score/affinity column. If None, auto-detected.
    label_column : str, optional
        Name of active/decoy label column. If None, auto-detected.

    Returns
    -------
    labels : np.ndarray
        Binary 1/0 array of active status.
    scores : np.ndarray
        Array of docking scores.
    properties_dict : dict
        Mapping property_name -> array of values.
    metadata : dict
        Dataset summary metadata.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    try:
        df = pd.read_csv(filepath, sep=None, engine='python', comment='#')
    except Exception:
        df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
        
    df.columns = [str(c).strip() for c in df.columns]
    
    # 1. Identify score column
    score_candidates = [
        "docking_score", "dock_score", "vina_score", "score", "affinity", "energy",
        "glide_gscore", "chemgauss4", "gold_fitness", "binding_affinity", "delta_g"
    ]
    if score_column is None:
        for cand in score_candidates:
            for col in df.columns:
                if cand == col.lower() or cand in col.lower():
                    score_column = col
                    break
            if score_column:
                break
                
    if score_column is None:
        # Fallback: pick the first numeric column with negative values or floating point
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                score_column = col
                break
                
    if score_column is None:
        raise ValueError(f"Could not identify a numeric score column in {filepath}")
        
    # 2. Identify label column (1/0, active/decoy, True/False)
    label_candidates = [
        "active", "label", "is_active", "class", "target", "activity", "category"
    ]
    if label_column is None:
        for cand in label_candidates:
            for col in df.columns:
                if cand == col.lower() or cand in col.lower():
                    label_column = col
                    break
            if label_column:
                break
                
    if label_column is None:
        # Check for binary numeric or string column
        for col in df.columns:
            if col != score_column:
                unique_vals = set(df[col].dropna().unique())
                if len(unique_vals) == 2 or unique_vals.issubset({0, 1, "0", "1", "active", "decoy", "Active", "Decoy", True, False}):
                    label_column = col
                    break

    if label_column is None:
        raise ValueError(f"Could not identify active/decoy label column in {filepath}")

    # Standardize labels to 1 (active) and 0 (decoy)
    raw_labels = df[label_column].to_numpy()
    binary_labels = np.zeros(len(raw_labels), dtype=int)
    for i, val in enumerate(raw_labels):
        val_str = str(val).strip().lower()
        if val_str in ["1", "1.0", "active", "true", "yes", "pos", "positive"]:
            binary_labels[i] = 1
        elif val_str in ["0", "0.0", "decoy", "inactive", "false", "no", "neg", "negative"]:
            binary_labels[i] = 0
        else:
            try:
                binary_labels[i] = 1 if float(val) > 0.5 else 0
            except ValueError:
                binary_labels[i] = 0
                
    scores = pd.to_numeric(df[score_column], errors='coerce').to_numpy(dtype=float)
    
    # 3. Extract physicochemical properties if present (MW, LogP, etc.)
    props_dict = {}
    known_props = ["mw", "molwt", "logp", "hbd", "hba", "rotbonds", "psa", "tpsa"]
    for col in df.columns:
        col_lower = col.lower()
        for kp in known_props:
            if kp in col_lower and col not in [score_column, label_column]:
                numeric_vals = pd.to_numeric(df[col], errors='coerce').to_numpy(dtype=float)
                props_dict[col] = numeric_vals
                break
                
    metadata = {
        "filename": os.path.basename(filepath),
        "score_column": score_column,
        "label_column": label_column,
        "n_total": len(scores),
        "n_actives": int(np.sum(binary_labels)),
        "n_decoys": int(len(binary_labels) - np.sum(binary_labels))
    }
    
    return binary_labels, scores, props_dict, metadata
