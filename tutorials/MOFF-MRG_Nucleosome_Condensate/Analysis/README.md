# Chromatin Condensate Trajectory Analysis Suite

This directory contains production-ready post-processing scripts designed to analyze 100-nucleosome coarse-grained chromatin slab trajectories. These tools measure the thermodynamic, structural, and network-level signatures driving liquid-liquid phase separation (LLPS).

---

## 1. Directory Blueprint
```text
analysis/
├── local_energy_rerun.py   # Quantifies individual nucleosome potential energy contributions
├── contacts_rerun.py       # Tracks structural inter-nucleosomal coordination profiles
└── network_rerun.py        # Generates frame-by-frame adjacency network contact maps

```

---

## 2. Script Architecture & Technical Details

### 🔋 `local_energy_rerun.py`

* **Core Physics:** This script isolates individual nucleosomes to evaluate their electrostatic screening footprint in the dense phase. It re-parses the targeted nucleosome's core and DNA atoms, dynamically remapping OpenMM atom names (`CA` $\rightarrow$ `NP` and `DN` $\rightarrow$ `ND`) to isolate specific force calculation groups.
* **Integrator Setup:** Computes a virtual potential energy rerun using a `NoseHooverIntegrator` via the `OpenABC` platform back-end.
* **Output:** Generates `openmm_energies.csv`, tracing the frame-by-frame electrostatic potential energy matrix ($100 \times N_{\text{frames}}$) in $\text{kJ/mol}$ using a custom Debye-Hückel switching map (`force_group=9`).

### 🧬 `contacts_rerun.py`

* **Core Physics:** Measures structural crowding by calculating coordination numbers between flexible histone tails and exposed terminal/linker DNA segments across the multi-array condensate.
* **Screening Parameters:** Calculates an explicit interaction cutoff distance calibrated directly to the ionic environment's Debye screening length ($\kappa$):

$$\text{Cutoff} = \kappa \times 10 = 0.6821\text{ nm} \times 10 = 6.821\text{ Å}$$


* **Optimization:** Leverages `MDAnalysis` to gather index offsets across all 100 nucleosomes, executing a broadcasted, pairwise distance matrix check across each trajectory frame while ignoring self-interaction noise.
* **Output:** Generates `contacts_200mM.csv` containing raw coordination metrics.

### 🕸️ `network_rerun.py`

* **Core Physics:** Treats the chromatin condensate as a topological network to quantify inter-nucleosomal connectivity distributions. Rows represent donor histone tail sets and columns represent acceptor DNA targets.
* **Graph Adjacency:** Evaluates structural connections between every unique pair combinations of nucleosomes $r$ and $c$ (where $r \neq c$).
* **Output:** Spits out a frame-by-frame adjacency matrix sequence (`contact_map_[frame]_frame.csv`). These files map the evolving interaction intensity across time, providing the baseline dataset used to construct the network robustness charts (**Figure S9**).

---

## 3. Configuration & Execution

Before initiating a run, ensure that your core path directory contains your primary structure and trajectory files (`cg_n_nucl.pdb` and `100_nucleosome_condensate_200mM.xtc`).

### Adjusting User Inputs:

For alternate environmental conditions, update the top parameter blocks in the files:

```python
# Set environment baseline (Example for 200mM)
salt = 200
kappa = 0.6821065688598205  # Corresponding Debye screening radius (nm)

```

### Execution Command:

```bash
python local_energy_rerun.py
python contacts_rerun.py
python network_rerun.py

```

---

## 4. Required Environments & Dependencies

Ensure your active Python environment features the following verified analytical packages:

* `openmm` (or `simtk.openmm`)
* `openabc`
* `mdtraj`
* `MDAnalysis`
* `numpy` & `pandas`

```
