<h1 align="center">Modelling mechanochemical coupling  in optogenetically activated cell layers</h1>

<h3 align="center">-- Finite Element Code -- </h1>

Implementation of the Finite Element Scheme for the publication 

Wörthmüller, D., Ziebert, F. & Schwarz, U. S. **Modeling mechanochemical coupling in optogenetically activated cell layers**. Biophysical Journal (2025)  [10.1016/j.bpj.2025.10.002](https://doi.org/10.1016/j.bpj.2025.10.002)

The simulations used to obtain the results in the publication are based on the open source finite element package [FEniCS](https://fenicsproject.org) and are implemented in python. 

The Dockerfile can be used to setup a docker image and container to get the latest stable version of legacy FEniCS. For instructions please follow the official [Docker instructions](https://docs.docker.com). 
The finite element meshes were created with the open source mesh generator [gmsh](https://gmsh.info). The final mesh files are included. 
Simulation output data have been excluded from the repository due to file sizes. Further, we have only included one single run for each condition. All other runs are identical up to the simulation parameters. 

Some of the data has been generated on run time and the remaining data has been obtained post simulation.

The yaml file is a template for a conda environment used for the analysis and plotting.
It does not contain any specific packages and can also be created from scratch.
To create the envirnoment run

```bash
conda env create -f FILENAME.yml
```

In the following we shortly summarize the main files. More details of the implementation can be found in the supplemental text of the publication. 

#### Doublet:

```text
Doublet/
├── create_csv.py                       #python file creating runParams.csv which contains paramater combinations for all runs
├── runParams.csv                       #parameter csv for all runs
├── run17/                              #one specific run (one line of runParams.csv)
│   │
│   ├── doublet_optogenetics.py         #main simulation script
│   ├── inputParams.csv                 #parameter file for main simulation script
│   ├── boundaries/subdomains.xdmf      #shows tags for cell-cell boundaries and cells
│   ├── analysis.py                     #single run analysis script to be executed after simulation is finished
│   ├── ...                             #mesh files .geo .xml .msh...
│   ├── Data/                           #containes data produced by doublet_optogenetics.py and analysis.py
│   └── Results/                        #contains plot directory produced by analysis.py
```

#### Cell chain:

```text
CellChain/
├── create_csv.py                       #python file creating runParams.csv which contains paramater combinations for all runs
├── runParams.csv                       #parameter csv for all runs
├── run0143/                            #one specific run (one line of runParams.csv)
│   │
│   ├── cell_chain_optogenetics.py      #main simulation script
│   ├── inputParams.csv                 #parameter file for main simulation script
│   ├── boundaries/subdomains.xdmf      #shows tags for cell-cell boundaries and cells
│   ├── kymograph_cell_line_multiple.py #single run analysis script to be executed after simulation is finished, creates kymograph plots
│   ├── ...                             # other files only relevant for multi-run analysis and mesh files .geo .xml .msh...
│   └── Data/                           #contains data produced by cell_chain_optogenetics.py and kymograph_cell_line_multiple.py
```

#### Tissue:

```text
Tissue/
├── tissue.py                           #main simulation script
├── inputParams.csv                     #parameter file for main simulation script
├── boundaries/subdomains.xdmf          #shows tags for cell-cell boundaries and cells
├── ...                                 #mesh files .geo .xml .msh...
└── Data/                               #contains data produced by tissue.py
```
