"""
Generates synthetic DUD-E-like screening CSV for testing and demonstration.
"""

import os
import numpy as np
import pandas as pd


def generate_screening_benchmark(output_path: str = "benchmark_screening.csv"):
    n_actives = 50
    n_decoys = 1500
    
    rng = np.random.default_rng(101)
    
    # Generate realistic compound IDs
    active_ids = [f"ACTIVE_{i+1:03d}" for i in range(n_actives)]
    decoy_ids = [f"DECOY_{i+1:04d}" for i in range(n_decoys)]
    
    # Docking affinities (kcal/mol)
    active_scores = rng.normal(-9.6, 0.9, size=n_actives)
    decoy_scores = rng.normal(-6.7, 1.1, size=n_decoys)
    
    # Molecular properties
    active_mw = rng.normal(390.0, 40.0, size=n_actives)
    decoy_mw = rng.normal(385.0, 45.0, size=n_decoys)
    
    active_logp = rng.normal(2.9, 0.6, size=n_actives)
    decoy_logp = rng.normal(2.8, 0.7, size=n_decoys)
    
    df_actives = pd.DataFrame({
        "Compound_ID": active_ids,
        "Docking_Score_kcal_mol": active_scores,
        "Activity_Status": ["active"] * n_actives,
        "Molecular_Weight": active_mw,
        "LogP": active_logp
    })
    
    df_decoys = pd.DataFrame({
        "Compound_ID": decoy_ids,
        "Docking_Score_kcal_mol": decoy_scores,
        "Activity_Status": ["decoy"] * n_decoys,
        "Molecular_Weight": decoy_mw,
        "LogP": decoy_logp
    })
    
    df_total = pd.concat([df_actives, df_decoys], ignore_index=True)
    # Shuffle
    df_total = df_total.sample(frac=1.0, random_state=42).reset_index(drop=True)
    
    df_total.to_csv(output_path, index=False)
    print(f"Benchmark screening dataset generated at: {os.path.abspath(output_path)}")


if __name__ == "__main__":
    generate_screening_benchmark()
