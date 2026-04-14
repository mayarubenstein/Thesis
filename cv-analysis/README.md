# Cyclic Voltammetry Analysis Tools

A Python library for parsing, analyzing, and baseline-correcting cyclic voltammograms (CVs), developed as part of my senior thesis in the Bocarsly Lab at Princeton University (2024).

> **Note:** Written as undergraduate research code in 2024. Uploaded as an archive of work I'm proud of.

## Background

Cyclic voltammograms are the primary experimental tool for studying electrochemical reactions. This library was built to automate the analysis of CV data from a study of chromium-gallium oxide electrodes for CO₂ reduction — specifically to enable comparison of experimental CVs against simulation output.

## What's in here

### `CVAnalysisFunctions.py`
A library built around a `CV` class that stores and analyzes CV data:

- **`importCV` / `importCVtoArray`** — parses raw CV output files into a numpy array of `CV` objects, handling variable electrodes and scan rates
- **`CV` class** — stores experimental conditions (electrode, analyte, scan rate) and measured data; methods include:
  - `Ehalf()` — estimates standard potential E½ by averaging peak potentials Ep1 and Ep2
  - `graph()` — plots the CV with labeled axes onto a given matplotlib subplot
  - `represent()` — generates a string description for titles/legends
  - `correctedEp1/2()`, `correctedMinCurrent/MaxCurrent()` — baseline-corrected peak finding (see below)

### Baseline Correction Algorithm
CVs have a sloped baseline from the electric double layer that obscures the true peak current. The correction algorithm:

1. Normalizes both current and potential arrays so the visual scale is consistent
2. Computes the discrete derivative dI/dE
3. Smooths with a Savitzky-Golay filter
4. Searches for the inflection point where |dI/dE| crosses 1 (transition from mostly-horizontal to mostly-vertical motion)
5. Uses the inflection point and endpoint to construct a linear baseline, then subtracts it

Validated against hand-analyzed data: Ep values within 5% error across all scan rates and electrodes, ip values within 5% with two outliers (one attributed to human error in the reference values).

### `testingBaseline.ipynb`
Demonstrates the library on ferricyanide [Fe(CN)₆]⁻⁴/⁻³ CVs across 3 electrodes (glassy carbon, two Cr-Ga oxide electrodes) at 9 scan rates (25–2000 mV/s). Includes Randles-Ševčík analysis showing diffusion-controlled behavior for the oxide electrodes.

## Dependencies

```
numpy
matplotlib
scipy  # for Savitzky-Golay filter
```

## Context

Senior thesis, Department of Chemistry, Princeton University, 2024. Bocarsly Lab.
