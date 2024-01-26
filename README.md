# Primary motor cortex (M1) circuits model
## Description
Multiscale model of mouse primary motor cortex (M1) developed using NetPyNE (www.netpyne.org).

The model is described in the following paper: 

- Dura-Bernal S, Neymotin SA, Suter BA, Dacre J, Moreira JVS, Urdapilleta E, Schiemann J, Duguid I, Shepherd GMG, Lytton WW. **Multiscale model of primary motor cortex circuits predicts in vivo cell type-specific, behavioral state-dependent dynamics.** Cell Reports (In Press) 


## Setup and execution

Requires NEURON with Python and MPI support. 

1. From /sim run `nrnivmodl ../mod`. This should create a directory called x86_64. 
2. To run type: `./runsim [num_proc]' or the equivalent `mpiexec -np [num_proc] nrniv -python -mpi init.py`

## Overview of file structure:

* /sim/init.py: Main executable; calls functions from other modules. Sets what parameter file to use.

* /sim/netParams.py: Network parameters

* /sim/cfg.py: Simulation configuration

* /cells: source .py, .hoc, .json or .pkl files for the different cell types used in the model; these will be imported into netpyne

* /conns: source .py, .json or .pkl files for the network connectivity; these will be imported into netpyne

* /mod: NMODL files containing the ionic channel and synaptic mechanisms used in the model 

* /data: where the model and simulation data is stored 


## Web-based graphical setup and execution via Open Source Brain (OSB) and NetPyNE GUI

1. Go to the OSB M1 model repository: https://v2.opensourcebrain.org/repositories/60

2. Log in to OSB with your username and password (you can easily register if you don't yet have an account)

3. Once inside the M1 repository, click on "Create New Workspace" and go to this new Workspace when prompted.

4. Once inside the Workspace, click on "Open with ... NetPyNE" (select option from drop-down menu)

5. Once inside the NetPyNE GUI, click on File --> Import --> From Python ...

6. Select the following options:
    - NetParams file: /home/jovyan/netpyne/workspace/Multiscale model of primary motor cortex (M1) circuits developed in NetPyNE/main/sim/netParams.py
    - NetParams variable: netParams
    - SimConfig file: /home/jovyan/netpyne/workspace/Multiscale model of primary motor cortex (M1) circuits developed in NetPyNE/main/sim/netParams.py
    - SimConfig variable: cfg
    - Path to mod files: /home/jovyan/netpyne/workspace/Multiscale model of primary motor cortex (M1) circuits developed in NetPyNE/main/sim/mod
    - Compile mod files: True

    Note: This step will be simplified in an upcoming version of OSB and NetPyNE GUI (early 2023), in which the model will be automatically loaded.

7. You can now graphically explore, simulate and analyze the model structure and outputs.
    
    Note: This model requires approximately 2 hours on 96 cores for 1 second of simulation. An upcoming version of OSB / NetPyNE GUI will soon allow you to submit simulations to the Neuroscience Gateway (NSG) and other HPCs. In this version, you can still inspect all the components of the model (populations, cells, connectivity etc) and run smaller scale simulations (e.g. by reducing the cell density of populations) within the up to 4 cores provided currently by OSB. 


For further information please contact: salvadordura@gmail.com 

