<h1 align="center">Curvature instability of an active gel growing on a wavy membrane</h1>

<h3 align="center">-- Finite Element Code -- </h1>

Implementation of the Finite Element Scheme for the publication 
K. Mihali, D. Wörthmüller, and P. Sens, **Curvature instability of an active gel growing on a wavymembrane**, arXiv:2510.17701 (2025) https://doi.org/10.48550/arXiv.2510.17701

The simulations used to obtain the results in the publication are based on the open source finite element package [FEniCS](https://fenicsproject.org) and are implemented in python. 

The Dockerfile can be used to setup a docker image and container to get the latest stable version of legacy FEniCS. For instructions please follow the official [Docker instructions](https://docs.docker.com). 
The finite element meshes are created with the python interface of the open source mesh generator [gmsh](https://gmsh.info). 
Simulation output data have been excluded from the repository due to file sizes. Further, we have only included one single run. All other runs are identical up to the simulation parameters. 

Some of the data is generated on run time and the remaining data is obtained post simulation.

The yaml file is a template for a conda environment used for the analysis and plotting.
It does not contain any specific packages and can also be created from scratch.
To create the envirnoment run

```bash
conda env create -f FILENAME.yml
```

In the following we shortly summarize the main files. More details of the implementation can be found in the supplemental text of the publication. 

```text
fem/
├── fem_main.py                       #python file to solve the PDE problem
├── inputParams.csv                   #parameter csv
├── Data/                             #simulation output -> xdmf, pvd, vtu and npz files
├── Meshfiles/                        #gmsh-files and other mesh-related xdmf files highlighting boundary tags
├── inputParams_add.csv               #a csv file with additionally on run-time computed values such as forces
├── analysis.py                       #python file to determine stress values on membrane surface -> written to inputParams_add.csv, visual checkups to ensure global force balance 
├── errMeasure.txt                    #error-measure value to inspect convergence 
├── pdf's and png's                   #graphical output from analysis.py for sanity checks
```
