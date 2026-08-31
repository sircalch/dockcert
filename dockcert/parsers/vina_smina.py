"""
Parsers for AutoDock Vina, Smina, and GNINA output logs.
"""

from typing import List, Dict, Any, Tuple
import os


def parse_vina_log(filepath: str) -> Tuple[List[float], List[float], List[float]]:
    """
    Parses AutoDock Vina / Smina output log files to extract binding affinity and RMSD values.

    Parameters
    ----------
    filepath : str
        Path to Vina log file.

    Returns
    -------
    affinities : list of float
        Estimated binding affinity (kcal/mol) for each mode.
    rmsd_lb : list of float
        RMSD lower bound from best mode.
    rmsd_ub : list of float
        RMSD upper bound from best mode.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    affinities = []
    rmsd_lb = []
    rmsd_ub = []
    
    in_table = False
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_s = line.strip()
            if line_s.startswith("-----+") or line_s.startswith("---+"):
                in_table = True
                continue
                
            if in_table:
                parts = [p.strip() for p in line_s.replace("|", " ").split()]
                if len(parts) >= 4 and parts[0].isdigit():
                    try:
                        aff = float(parts[1])
                        lb = float(parts[2])
                        ub = float(parts[3])
                        affinities.append(aff)
                        rmsd_lb.append(lb)
                        rmsd_ub.append(ub)
                    except ValueError:
                        break
                elif len(parts) > 0 and not parts[0].isdigit() and len(affinities) > 0:
                    break
                        
    return affinities, rmsd_lb, rmsd_ub
