```markdown
# Chromatin Condensate Slab Simulation Engine

This folder contains the core configuration and execution scripts used to build and simulate 100-nucleosome chromatin arrays within an elongated slab geometry. The setup utilizes the `OpenABC` platform to implement the coarse-grained **MOFF** protein potential and **MRG** double-stranded DNA model.

---

## 1. Directory Structure

Ensure your local folder maintains the following architecture before executing the pipeline:

```text
simulation/
├── Setup.py                # Main workflow script (System construction & NVT run)
├── nucl_utils.py           # Core library containing rigid body mapping utilities
└── starting_structure/
    ├── histone.pdb         # Atomistic/Coarse-grained reference structural core
    └── dna.pdb             # Atomistic/Coarse-grained baseline tracking DNA coordinates

```

---

## 2. Technical & Physical Specifications

* **Force Field Framework:** Operates at a unified residue-resolution level. Histones are treated using the C$\alpha$-level **MOFF** model, while DNA parameters employ the **MRG** model.
* **Calibrated DNA Rigidity:** Linker and wrapped DNA features a tuned `bonded_energy_scale = 0.8`. This mechanical calibration is explicitly validated to reproduce the empirical **50 nm persistence length** of naked double-stranded DNA.
* **Slab Phase Coexistence Geometry:** Molecules are packed into an anisotropic, elongated box container ($49 \times 49 \times 300\text{ nm}$). The extended $z$-axis provides the thermodynamic space needed for a dense chromatin core to phase-separate and continuously coexist alongside a surrounding dilute phase.
* **Rigid Body Constraints:** To maintain structural integrity during dense crowding conditions, each nucleosome enforces a rigid body grouping consisting of the central histone octamer and the internal 73 base pairs of wrapped DNA.

---

## 3. Configuration Parameters

The baseline environment variables are hardcoded within `Setup.py` to target physiological conditions:

* **Salt Environment:** $200\text{ mM}$ monovalent salt concentration (governed by custom Debye-Hückel switching map electrostatics).
* **Thermal Controls:** $300\text{ K}$ regulated via a `LangevinMiddleIntegrator`.
* **Kinetics Acceleration:** The integrator friction coefficient is set to a low value ($0.01\text{ ps}^{-1}$) with a **10 fs timestep** to accelerate conformation sampling.
* **Hardware:** Default computational target is configured for **CUDA** GPU acceleration.

---

## 4. Pipeline Execution & Outputs

### Running the Workflow

To initiate the nucleosome packing, structural assembly, force group assignments, and subsequent 50 ns ($5,000,000$ steps) production NVT compression trajectory, execute:

```bash
python Setup.py

```

### Generated Outputs

Upon successful completion, the script will populate your workspace with four key runtime files required for downstream post-processing analysis:

1. `cg_n_nucl.pdb`: The full assembled topology file holding the initial coordinates of all 100 packed nucleosomes.
2. `system.xml`: Serialized system parameter layout holding active force groups and constraints.
3. `traj_temp300.dcd`: Binary coordinate trajectory file recording periodic matrix frames every 50,000 integration steps.
4. `state_report.txt`: Tab-delimited performance report logging computational steps, absolute temperatures, potential energy trends, and processing velocities.
5. `state_cpt.xml`: Final thermodynamic checkpoint holding coordinate and velocity vectors, utilized to resume trajectories smoothly without state loss.

```

```
