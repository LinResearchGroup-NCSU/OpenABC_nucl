# Chromatin Condensate Trajectory Analysis Suite

This directory contains production-ready post-processing scripts designed to analyze 100-nucleosome coarse-grained chromatin slab trajectories. These tools measure the thermodynamic, structural, and network-level signatures driving liquid-liquid phase separation (LLPS).

---

## 1. Directory Blueprint
```text
analysis/
├── local_energy_rerun.py   # Quantifies individual nucleosome potential energy contributions
├── contacts_rerun.py       # Tracks structural inter-nucleosomal coordination profiles
└── network_rerun.py        # Generates frame-by-frame adjacency network contact maps
