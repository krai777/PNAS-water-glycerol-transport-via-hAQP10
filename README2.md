#Monomer_permeation
A Python-based analysis tool for detecting and quantifying molecular permeation events through membrane protein pores using MD trajectories analyzed with MDAnalysis.
This script tracks molecules (e.g. water or glycerol) crossing a cylindrical pore region relative to the center of mass (COM) of individual protein chains.
#Features
*.Detects full permeation events through a pore
*.Tracks permeation directionality:
  lower → upper (up)
  upper → lower (down)
*.Supports multiple independent protein chains
*.Uses a cylindrical geometric criterion
*.Stores:
  frame-by-frame permeation counts
  cumulative permeation statistics
  individual permeation event details
  Compatible with wrapped MD trajectories
