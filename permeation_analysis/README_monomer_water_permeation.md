# Monomer_Permeation_Tracker

A Python tool for tracking and counting molecular permeation events through a membrane protein pores from MD simulation trajectories.

## Overview

`Monomer_Permeation_Tracker` analyzes molecular dynamics (MD) trajectories to detect when water molecules (or other permeating species like glycerol) pass through a cylindrical pore region of a membrane protein. It tracks each permeation event per monomer chain, records directionality (up/down), and outputs per-frame counts and per-event details.

Designed for use with homotetrameric channel proteins (e.g. aquaporins), with support for four monomer chains (A–D).

## Requirements

- Python 3.7+
- [MDAnalysis]([https://www.mdanalysis.org/])
- NumPy

Install dependencies:

```bash
pip install MDAnalysis numpy
```

## Usage

Edit the `__main__` block at the bottom of the script with your file paths and parameters:

```python
psf_file = "path/to/your.psf"
dcd_file = "path/to/your.dcd"

tracker = PermeationTracker(
    psf_file=psf_file,
    dcd_file=dcd_file,
    cylinder_radius=8.0,   # Å, radial cutoff for pore region
    z_min=-18.0,           # Å, lower z boundary relative to chain COM
    z_max=18.0             # Å, upper z boundary relative to chain COM
)

tracker.run_analysis()
```

Then run:

```bash
python permeation_tracker.py
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `cylinder_radius` |  Radial cutoff defining the cylindrical pore |
| `z_min` |  Lower z boundary relative to chain center of mass |
| `z_max` |  Upper z boundary relative to chain center of mass |

## Output Files

For each chain (A, B, C, D), two files are written:

**`NEWwatcounts_{chain}.txt`** — per-frame permeation counts:
```
FRAME   up_to_down   down_to_up   total_permeation   cumsum_permeation
```

**`NEWwatdetails_{chain}.txt`** — per-event details:
```
MoleculeID   entry_frame   exit_frame
```

## Permeating Species

By default the script tracks TIP3 water oxygens. To track a different species (e.g. glycerol), update the selection in `__init__`:

```python
# Water (default)
self.waters = self.u.select_atoms("resname TIP3 and name OH2")

# Glycerol example
self.waters = self.u.select_atoms("resname MGLY and name C1")
```

## Notes

- A permeation event is counted only when a molecule enters one side of the pore, crosses the central plane (z = 0 relative to COM), and exits from the opposite side.
- Each chain's pore axis is defined relative to that chain's center of mass, making the analysis robust to global translations.
