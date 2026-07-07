# Coarse-Grained Chromatin Slab Simulation (OpenABC MOFF/MRG Forcefield)

This repository provides an automated workflow to set up, parameterize, and execute high-throughput molecular dynamics trajectories of multi-nucleosome arrays inside an elongated slab geometry. The software is optimized to investigate liquid-liquid phase separation (LLPS), network connectivity, and thermodynamic binodal boundaries of chromatin as a function of ionic strength, DNA mechanics, and epigenetic tailoring.

## 1. Physical & Architectural Overview
* **Forcefield Framework:** Employs a unified residue-resolution model combining the **MOFF** (Molecular model for Open chromatin with Fine-tuned Flexible histone tails) protein potential with the **MRG** (Multi-Scale Rigid-Body/Gaussian) double-stranded DNA model via the `OpenABC` API.
* **Mechanical Validation:** The model utilizes a calibrated DNA bending rigidity parameter (`bonded_energy_scale=0.8`) which accurately reproduces the empirical **50 nm persistence length** characteristic of standard 150-bp naked DNA segments.
* **System Geometry:** Simulations are housed within an anisotropic slab container ($49 \times 49 \times 300\text{ nm}$). The extended $z$-axis allows a stable dense chromatin condensate core to spontaneously form and remain in continuous coexistence with a surrounding dilute phase.
* **Rigid Body Constraints:** To maintain structural realism, each individual nucleosome features a rigid core containing the histone octamer bound alongside the internal 73 base pairs of wrapped DNA.

## 2. Directory & Input Requirements
Before running the primary simulation script, verify that your working environment contains the following file structure:
```text
├── run_chromatin_slab.py       # Main setup and integration script
├── nucl_utils.py               # Auxiliary chromatin topology generator functions
├── starting_structure/
│   ├── histone.pdb             # Atomistic/Coarse-grained baseline structure of single histone core
│   └── dna.pdb                 # Atomistic/Coarse-grained base map of single nucleosome DNA track
└── equil-200mM/
    └── output_NPT.dcd          # Pre-compressed equilibration coordinate track
