# Nanopore Electrode Simulation

A modular Python simulation framework for modelling electrochemical reactions at porous insulating electrodes, developed as my senior thesis in the Bocarsly Lab at Princeton University (2024). The scientific motivation is understanding the surprising catalytic activity of chromium-gallium oxide nanopores for CO₂ reduction.

> **Note:** This is research code from 2024. The core data structures (`volumeHT`, `boxes`) are fully implemented and tested. The electrode and simulation modules are implemented but not fully debugged — this was ongoing work at the time of submission. Uploaded as an honest archive of the algorithmic design work.

## Motivation

Insulating nanoporous electrodes (Cr₂O₃/Ga₂O₃) show anomalously high faradaic efficiency for CO₂ reduction. The hypothesis is that mechanical confinement — not chemical or quantum effects — drives this. This simulation framework is designed to test that by modeling particle diffusion and reaction in explicitly-specified pore geometries.

## Architecture

The simulation is built around five loosely-coupled modules:

```
volumeHT          ← 3D spatial hash table for collision checking
    └── boxes     ← Box variants (radiusBox, overlapBox, augmented versions)
electrode         ← Flat and porous electrode modules
    └── pore      ← Cylindrical pore geometry, collision/rebound logic
molecule          ← Chemical species with diffusion and reactivity
simulation        ← Outer loop, integrates all modules
```

### `boxes.py` — Box variants
Four box types with a shared API (`add`, `delete`, `population`, `particles`, `allowedMoveInBox`, `allowedMoveOutOfBox`), allowing the client to swap them interchangeably in `volumeHT`:

- **`radiusBox`** — hard-sphere collision checking using `SortedList` for O(log n) deletion
- **`overlapBox`** — particles may overlap; used for the bulk region where collision checking is unnecessary
- **`augRadiusBox` / `augOverlapBox`** — augmented variants that additionally track per-species particle counts

### `redoxSimulationCopy.py` — `volumeHT` and core simulation
**`volumeHT`** is the algorithmic centerpiece. It partitions 3D space into a dictionary of boxes for efficient collision checking:

- Computes integer box counts per dimension from a target particle density
- `coordsToBox()` — maps a coordinate to its box in O(1)
- `adjacent()` — returns the up-to-26 neighboring boxes of a given box
- `put()` — inserts a particle with collision checking; uses `fromCenter()` to skip adjacency checks when the particle is far enough from box edges (geometric optimization that avoids unnecessary neighbor traversal)
- `attemptMove()` — moves a particle with rehashing and automatic rollback on collision

The `molecule` class models [Fe(CN)₆]⁻³/⁻⁴ with Butler-Volmer reaction kinetics, diffusion constants, and hard-sphere radii. The `reactorTracker` class logs reaction events per timestep for CV reconstruction.

### `electrode.py`
Two electrode modules with a shared interface (designed to be swapped into the simulation with minimal code changes):

- **`electrode2D`** — flat conductive electrode; handles particle interception using parametric line equations, rebound via polar coordinates, and surface sticking via a Bernoulli parameter
- **`porousElectrode`** — porous electrode; delegates to `pore` and `areaHT` objects. `areaHT` is a 2D analogue of `volumeHT` that maps xy-boxes to the pores they intersect, for efficient pore-entry checking

The `pore` object implements cylindrical geometry with separate `floorReact`/`wallReact` booleans (insulating walls, conducting floor), wall-collision detection via quadratic formula on the parametric trajectory, and wall-rebound via a rotation matrix derived from the tangent plane at the particle's position.

## What works / what doesn't

| Component | Status |
|---|---|
| `volumeHT` + `boxes` | ✅ Fully implemented and tested |
| `molecule`, `reactorTracker` | ✅ Implemented |
| `electrode2D`, `porousElectrode` | ⚠️ Implemented, not fully debugged |
| Full CV simulation | ⚠️ Not completed |

## Dependencies

```
numpy
sortedcontainers
matplotlib
```

## Context

Senior thesis, Department of Chemistry, Princeton University, 2024. Bocarsly Lab (PI: Prof. Andrew Bocarsly). Research question: can the catalytic activity of insulating Cr-Ga oxide nanopores for CO₂ reduction be explained purely by mechanical confinement effects?
