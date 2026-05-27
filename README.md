# water-permeation-analysis
This repository contains a VMD/Tcl script for counting water, ion, and small molecules permeation events via membrane/protein channel from molecular dynamics unwrapped trajectories.
The script analyzes a trajectory using a PSF and DCD file and writes the frame-wise permeation counts to `watnum.dat`.
## Purpose

The script is designed to quantify molecular permeation across a membrane/channel system by monitoring molecules that cross from one side of the simulation box to the other. It currently counts:

- water molecules (`OH2` atoms)
- sodium ions (`SOD`)
- chloride ions (`CLA`)
- glycerol molecules (`MGLY`, atom name `C2`)
  
You can change according to your system specifications.
The output can be used to calculate net permeation events, cumulative transport, for molecular dynamics analysis.

## Preprocessing requirement
Before using this script, the input structure must be prepared with appropriate segment names that identify molecules initially located on the upper leaflet side and lower leaflet side of the system.
In particular, the PDB/PSF files must contain segment names such as:

- `WATUP` and `WATDOWN` for water
- `SODUP` and `SODDOWN` for sodium ions
- `CLAUP` and `CLADOWN` for chloride ions
- `GLYUP` and `GLYDOWN` for glycerol

This is essential because the script detects permeation using selections. Therefore, the PSF and PDB files should be generated or modified before the analysis so that molecules on opposite sides of the membrane/channel are assigned the correct segment names.

## Usage
Run the script using VMD in text mode:

```bash
vmd -dispdev text -e water_permeation.tcl -args system.psf trajectory.dcd
```
where:

- `system.psf` is the PSF file
- `unwrap_trajectory.dcd` is the unwrapped trajectory file
- `water_permeation.tcl` is the analysis script

## Notes

- The script assumes that the membrane/channel axis is aligned with the z-axis.
- The crossing boundaries are currently defined using `z < -15` and `z > 15`.
- These values may need to be adjusted depending on the system geometry.
- The segment names in the PSF/PDB must match the names used in the script.
