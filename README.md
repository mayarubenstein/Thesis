# Atomistic Simulations for Insulating Nanoporous Electrodes

Senior thesis, Department of Chemistry, Princeton University, 2024.  
Bocarsly Lab (PI: Prof. Andrew Bocarsly).

> **Note:** This is undergraduate research code from 2024, uploaded as an honest archive. Some modules are fully working; others were works in progress at the time of submission. See each subfolder's README for details.

## Research Question

Chromium-gallium oxide (Cr₂O₃/Ga₂O₃) nanoporous electrodes show anomalously high faradaic efficiency for CO₂ reduction — despite being electrical insulators. The hypothesis investigated here is that **mechanical confinement**, not chemical or quantum effects, is responsible. This project builds the simulation infrastructure needed to test that hypothesis.

## Repository Structure

### [`cv-analysis/`](./cv-analysis/)
A working Python library for parsing, analyzing, and baseline-correcting cyclic voltammograms (CVs). Used to process experimental ferrocyanide CV data from the Bocarsly Lab and validated against hand-analyzed results.

### [`simulation/`](./simulation/)
A modular Python simulation framework for modeling electrochemical reactions at porous electrodes via particle diffusion. The centerpiece is `volumeHT`, a 3D spatial hash table for efficient collision checking. Also includes flat and porous electrode modules with cylindrical pore geometry.

## Dependencies

```
numpy
matplotlib
scipy
sortedcontainers
```
