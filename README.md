# Restriction decoder for the domain wall colour code

This repository contains code for code capacity simulations of the X3Z3 domain wall colour code with a restriction decoder. 
The implemented decoder is explained in section 2.2 and 2.3 of https://iopscience.iop.org/article/10.1088/1367-2630/ab68fd/meta.
I use PyMatching for performing MWPM. 

### Run Tests

The following assumes you've cloned the repo and your current working directory is the repo root.

```bash
# (in a fresh virtual environment at repo root)
pip install -r requirements.txt
cd src
python3 -m pytest
```

### Running simulations and generating plots

The file src/run_calculate_ler.py can be used to run simulations.
The data used to create the threshold plot in the paper is in the directory data_18_11
Plots are created in plots.create_threshold_plot.ipynb