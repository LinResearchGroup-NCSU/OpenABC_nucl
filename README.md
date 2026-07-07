# Chromatin Phase Separation Workflow (OpenABC Nucleosome Condensate Engine)

This repository provides an automated, high-throughput molecular dynamics workflow designed to simulate **100-nucleosome chromatin arrays** inside an elongated slab geometry. Built upon the `OpenABC` platform, this specific package integrates the **MOFF** protein potential and the **MRG** coarse-grained DNA model to investigate the thermodynamic, mechanical, and epigenetic mechanisms governing liquid-liquid phase separation (LLPS) of genome architecture.

---

## Highlighted Tutorial

> 📖 **Getting Started:** We provide a comprehensive, step-by-step pipeline demonstration in our **Nucleosome Condensate Tutorial**. This python file guides users through structural setup, custom force integrations, and production NVT trajectory executions. 
> 
> **Go straight to the pipeline:** `tutorials/MOFF-MRG_Nucleosome_Condensate/Simulation/Setup.py`

---

## Features

- **Slab Phase Coexistence:** Automated assembly of $49 \times 49 \times 300\text{ nm}$ simulation boxes tailored to capture clean dense/dilute phase boundary profiles.
- **Calibrated DNA Rigidity:** Implementation of a tuned DNA persistence model (`bonded_energy_scale=0.8`) matching the empirical 50 nm baseline length of native linker DNA.
- **Post-Processing Analytics:** Complete standalone suites for tracking localized electrostatic potential energy, structural coordination numbers, and graph-theory network adjacency profiles.

---

## Installation & Environment Setup

This implementation operates on top of the main **OpenABC** library. For detailed background information, comprehensive class manuals, and underlying code architectures, please refer directly to the official [OpenABC Repository](https://github.com/ZhangGroup-MITChemistry/OpenABC).

### 1. Environment Requirements
Create a virtual environment with `conda` (or `miniconda`) incorporating GPU-accelerated OpenMM and basic structural analysis libraries:
```bash
conda create -n chromatin_env -c conda-forge openmm mdtraj mdanalysis numpy pandas python=3.10
conda activate chromatin_env
```

### 2. Install the OpenABC Core API

Install the underlying forcefield parser package via pip:

```bash
pip install openabc

```

### 3. Deploy this Repository

Clone this branch to access the specific nucleosome slab assembly and analysis pipeline:

```bash
git clone [https://github.com/YourUsername/Your_Chromatin_Repository.git](https://github.com/YourUsername/Your_Chromatin_Repository.git)
cd Your_Chromatin_Repository

```

---

## Usage Workflow

The workspace is divided into a **Production Execution** phase and a subsequent **Trajectory Post-Processing** phase.

### Phase 1: Running the 100-Nucleosome Slab Simulation

The main production script initializes the anisotropic box volume, enforces rigid-body constraints across the 73 bp wrapped histone cores, maps the Debye-Hückel electrostatics screening map, and executes a 50 ns production NVT sampling run.

Ensure `histone.pdb` and `dna.pdb` are placed in your `starting_structure/` folder, then run:

```bash
python Setup.py

```

*Outputs generated:* `traj_temp300.dcd` (trajectory tracks) and `state_report.txt` (thermodynamic logs).

### Phase 2: Post-Processing & Validation Analytics

Navigate to the `analysis/` folder to process your simulation trajectories using the specialized scripts mapped out below:

1. **Local Energy Evaluation:** Remaps localized atom groups to evaluate virtual electrostatic potentials across the dense core interface.
```bash
python analysis/local_energy_rerun.py

```


2. **Inter-Nucleosomal Coordination:** Computes structural contact timelines between flexible tail segments and exposed DNA segments based on ionic Debye screening bounds.
```bash
python analysis/contacts_rerun.py

```


3. **Graph-Theory Network Maps:** Generates frame-by-frame adjacency matrices monitoring topological network robustness.
```bash
python analysis/network_rerun.py

```



*For granular details on input array dimensions and index configurations within these python codes, check the local `analysis/README.md` documentation.*

---

## References & Citations

If you utilize this workflow, forcefield architecture, or analysis suite in your research, please cite the primary OpenABC platform alongside the structural models:

* **OpenABC Framework:** Move, et al. "OpenABC Enables Flexible, Simplified, and Efficient GPU Accelerated Simulations of Biomolecular Condensates." *Bioinformatics*, 2023. DOI: [10.1101/2023.04.19.537533](https://doi.org/10.1101/2023.04.19.537533)
* **MOFF + MRG Model:** Latham, A. P., & Zhang, B. "On the stability and layered organization of protein-DNA condensates." *Biophysical Journal*, 121(9), 1727-1737, 2022. DOI: [10.1016/j.bpj.2022.03.024](https://www.google.com/search?q=https://doi.org/10.1016/j.bpj.2022.03.024)

```

```


