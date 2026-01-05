# CHAOS database

**CHAOS (Computed High-Accuracy Observables and Sigma-profiles)** is a large-scale, internally consistent quantum-chemical database providing **σ-profiles and a broad set of additional molecular descriptors for 53,091 molecules**.  
All quantum-chemical entries were generated with a standardized workflow at the **ωB97X-D/def2-TZVP** level of theory, ensuring strict comparability across the entire dataset.

---

The full CHAOS data release is hosted on **Zenodo**:

- **Zenodo DOI:** https://doi.org/10.5281/zenodo.17691925

This GitHub repository will provide documentation, and scripts soon:

- **GitHub:** https://github.com/GondLTD/CHAOS_database

**License:** CC-BY-4.0 — data may be used, modified, and redistributed with attribution.

## Data format

Each molecule is stored as a structured JSON object with six main blocks:

- `general`
- `structural`
- `electronic`
- `vibrational`
- `NMR`
- `solvation`

Below is an overview of all entries in the CHAOS database, including computed properties.  
For each property we list the JSON key followed by a short description.
There is also one JSON file for all molecules available on Zenodo, omitting the COSMO segment list for faster loading.

---

### General information

- `CanonicalSMILES` — canonical SMILES of the molecule, canonicalized with RDKit version `2024.3.1` in Python 3.12.4.
- `MolecularMass` — molecular mass read from Gaussian output, float scalar in atomic mass units amu.
- `AtomList` — list of all atoms in the molecule. Each entry is a dictionary with keys `index`, `element`, and `atomic_number`.
- `not_converged` — Boolean flag for geometry optimization convergence. If negative frequencies persist after the second DFT cycle, set to `true`, otherwise `false`.

---

### Structural properties

- `Coordinates` — Cartesian coordinates in Gaussian standard orientation. Nested list; outer list maps to atoms ordered as in `AtomList`. Inner lists are `[x, y, z]` in Å.
- `Coordniates_Input` — Cartesian coordinates in Gaussian input orientation. Same structure/units as `Coordinates`. Available only for molecules with <52 atoms (incl. H). If not obtainable, entry is `None`; then use `Coordinates`.
- `RotConst` — rotational constants from Gaussian, list `[A, B, C]` in GHz with spectroscopic convention `A ≥ B ≥ C`. For linear molecules: list with single entry `[B]`.
- `RotTemps` — rotational temperatures from Gaussian, list `[T_A, T_B, T_C]` in K with convention `T_A ≥ T_B ≥ T_C`. For linear molecules: list with single entry `[T_B]`.
- `PointGroup` — Schönflies point group symbol as string (Gaussian notation).
- `SymNumber` — external symmetry number from Gaussian, integer.

---

### Electronic properties

- `Charge` — total molecular charge from RDKit, integer in elementary charges.
- `Multiplicity` — spin multiplicity. Computed by RDKit assuming high-spin for unpaired electrons from `GetNumRadicalElectrons`. Multiplicity of oxygen set to 3.
- `DipoleMoment` — scalar dipole magnitude from Gaussian, float in D.
- `DipoleMoment_au` — scalar dipole magnitude from Gaussian in atomic units (e·a₀).
- `DipoleMoment_SI` — scalar dipole magnitude from Gaussian in C·m.
- `DipoleVec` — dipole vector from Gaussian (input orientation), list `[μx, μy, μz]` in D.
- `DipoleVec_au` — dipole vector from Gaussian (input orientation), list `[μx, μy, μz]` in e·a₀.
- `DipoleVec_SI` — dipole vector from Gaussian (input orientation), list `[μx, μy, μz]` in C·m.
- `QuadrupoleMoment` — scalar magnitude of quadrupole matrix from Gaussian in D·Å, computed from traceless quadrupole diagonal elements.
- `QuadrupoleMat` — quadrupole matrix from Gaussian, dictionary with keys `xx`, `yy`, `zz`, `xy`, `xz`, `yz` in D·Å.
- `QuadrupoleMat_traceless` — traceless quadrupole matrix from Gaussian, dictionary with keys `xx`, `yy`, `zz`, `xy`, `xz`, `yz` in D·Å.
- `PolarMat` — dipole polarizability matrix from Gaussian, dictionary with keys `xx`, `yy`, `zz`, `xy`, `xz`, `yz` in (a₀)³.
- `PolarMat_SI` — dipole polarizability matrix from Gaussian in C²·m²·J⁻¹.
- `PolarMat_esu` — dipole polarizability matrix from Gaussian in cm³ (esu units).
- `PolarIso` — isotropic polarizability average in (a₀)³.
- `PolarIso_SI` — isotropic polarizability average in C²·m²·J⁻¹.
- `PolarIso_esu` — isotropic polarizability average in cm³.
- `PolarAniso` — polarizability anisotropy in (a₀)³.
- `PolarAniso_SI` — polarizability anisotropy in C²·m²·J⁻¹.
- `PolarAniso_esu` — polarizability anisotropy in cm³.
- `PartChargeMulliken` — per-atom Mulliken partial charges from Gaussian, list in elementary charges.
- `PartChargeAPT` — per-atom APT partial charges from Gaussian, list in elementary charges.
- `PartChargeMullikenHeavy` — Mulliken partial charges for heavy atoms with attached H charges summed onto bonded heavy atoms, list in elementary charges.
- `PartChargeAPTHeavy` — APT partial charges for heavy atoms with attached H charges summed onto bonded heavy atoms, list in elementary charges.
- `SCFEnergy` — self-consistent field energy from Gaussian, float in Hartree (H). Reference state is free electrons and atomic cores at infinite separation.
- `HOMOEnergy` — HOMO orbital energy from Gaussian, float in H.
- `LUMOEnergy` — LUMO orbital energy from Gaussian, float in H.
- `HLG` — HOMO–LUMO gap (LUMO − HOMO), float in eV.
- `HLG_Hartree` — HOMO–LUMO gap, float in H.

---

### Vibrational properties

- `Frequencies` — harmonic frequencies from Gaussian, list of floats in cm⁻¹. Length is `3N − F` with `F=5` for linear and `F=6` otherwise.
- `IRIntensities` — IR intensities corresponding to `Frequencies`, list in km·mol⁻¹.
- `ForceConstants` — force constants corresponding to `Frequencies`, list in mDyne·mol⁻¹.
- `ReducedMasses` — reduced masses corresponding to `Frequencies`, list in u.
- `ZPE` — zero-point energy from Gaussian, float in J·mol⁻¹. For anharmonic correction, scale frequencies by 0.95461 (Kesharwani et al.).
- `HeatCap` — molar isochoric heat capacity at 298 K from Gaussian, float in cal·mol⁻¹·K⁻¹. Same scaling recommendation as for ZPE.
- `Entropy` — absolute entropy at 298 K from Gaussian, float in cal·mol⁻¹·K⁻¹. Same scaling recommendation as for ZPE.
- `E_Thermal` — thermal energy at 298 K from Gaussian, float in kcal·mol⁻¹. Same scaling recommendation as for ZPE.

---

### NMR properties

- `ShieldIso` — per-atom isotropic shielding constants from Gaussian, list in ppm.
- `ShieldAniso` — per-atom shielding anisotropies from Gaussian, list in ppm.
- `SuscDia` — diamagnetic susceptibility tensor from Gaussian, dictionary with keys `xx`, `yy`, `zz`, `xy`, `xz`, `yz` in α²·a₀³.
- `SuscPara` — paramagnetic susceptibility tensor from Gaussian, dictionary with keys `xx`, `yy`, `zz`, `xy`, `xz`, `yz` in α²·a₀³.
- `SuscDiaIso` — isotropic average of diamagnetic susceptibility tensor, float in α²·a₀³.
- `SuscParaIso` — isotropic average of paramagnetic susceptibility tensor, float in α²·a₀³.

---

### Solvation properties

- `CavArea` — cavity surface area from Gaussian, float in (a₀)².
- `CavVol` — cavity volume from Gaussian, float in (a₀)³.
- `CavSegments` — number of COSMO surface segments from Gaussian, integer.
- `SCFCOSMOEnergy` — SCF energy of the solute inside the cavity from Gaussian, float in\*H. Note: internal energy, not Gibbs energy.
- `DielectricCorr` — dielectric correction term from Gaussian, float in kcal·mol⁻¹.
- `CavEnergy` — cavitation energy contribution from Gaussian, float in kcal·mol⁻¹.
- `DispersionEnergy` — dispersion contribution from Gaussian, float in kcal·mol⁻¹.
- `PauliRepulsion` — Pauli repulsion contribution from Gaussian, float in kcal·mol⁻¹.
- `Sigma_NHB` — normalized σ-profile for non–H-bonding atoms, Å², COSMO-SAC-dsp protocol. Bins from −0.025 to 0.025 e/Å² in steps of 0.001 e/Å².
- `Sigma_OH` — normalized σ-profile for hydroxyl H-bonding atoms, Å², COSMO-SAC-dsp protocol. Same binning.
- `Sigma_OT` — normalized σ-profile for other H-bonding atoms (non-OH), Å², COSMO-SAC-dsp protocol. Same binning.
- `Sigma_total` — unnormalized total σ-profile, Å², COSMO-SAC-dsp protocol. Same binning.
- `Norm_Sigma_total` — integral of total σ-profile; equals `CavArea` in Å².
- `Fraction_PartSigmas` — proportions of partial σ-profiles to total, list ordered `[NHB, OH, OT]`.
- `SegmentList` — full COSMO surface segment output from Gaussian, nested list; each segment entry contains:  
  `[segment_number, atom_index, x, y, z, charge, area, sigma, potential]`.  
  Provided **only in molecule-specific JSON files**.
