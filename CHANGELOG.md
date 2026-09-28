# Changelog

## 1.1.0 (2026-09-27)

Validated against real AutoDock Vina 1.2.7 redocking of trypsin–benzamidine (3PTB), streptavidin–biotin
(1STP) and HIV-1 protease–indinavir (1HSG) (`validation/redock_vina.py`). The redocking RMSDs now match
RDKit `CalcRMS` to within 0.006 Å for all 27 poses. Version 1.0.0 was off by up to 2.3 Å.

### Fixed
- **Multi-pose output.** A Vina/Smina PDBQT file was read as a single molecule containing all MODELs
  (81 atoms for a 9-atom ligand). As a result, the RMSD compared the reference with an arbitrary subset
  of atoms. Poses are now split by MODEL/ENDMDL (PDB/PDBQT) and `$$$$` (SDF).
- **AutoDock atom types.** Types in the PDBQT type column (`A`, `OA`, `NA`, `SA`, `HD`) were taken
  as element symbols. They are now mapped to elements, and polar hydrogens (`HD`) are excluded from
  the heavy-atom RMSD.
- **Symmetry correction.** The per-element Hungarian assignment ignored bonding. It gave a lower bound,
  not the RMSD, and depended on the reference atom order. The RMSD is now minimised over element- and
  bond-preserving graph isomorphisms, which is the definition used by RDKit and spyrmsd. The Hungarian
  assignment is kept only as a fallback when the two graphs differ.
- **Redocking verdict.**
  - Single run: the status was set by the best of all poses, so a wrong top-ranked pose could PASS.
    The verdict now uses the top-scored pose. A lower-ranked pose within 2 Å gives a WARNING
    (scoring failure), and the best of N poses is reported separately.
  - Benchmark of complexes: one lucky target no longer gives a PASS. The verdict is based on the
    success rate (PASS ≥ 70 %, WARNING ≥ 50 %).
- The generated methods text says the protocol "was tested" by redocking, not "validated",
  and reports the top-1 RMSD.

- **Enrichment factor cutoff.** The top x % now contains ceil(x·N) compounds, as in RDKit's
  `CalcEnrichment`; before it was round(x·N). EF_max is computed with the same cutoff. The ROC-AUC,
  BEDROC, RIE and EF values match RDKit exactly on tie-free scores (`tests/test_enrichment_rdkit.py`).

### Added
- `load_poses()`, `infer_bonds()`, the `elements_dock` argument of
  `calculate_symmetry_corrected_rmsd`, the `ranked_poses` argument of `evaluate_redocking_success`,
  and `rmsd_ranked_poses` in `assess_docking_quality`.
- CLI: `--docked-pose` reads every pose and `--vina-log` lists the affinities next to each pose's RMSD.
- Tests on the real redocking outputs (`tests/test_redock_real.py`).
