"""
init.py

Starting script to run NetPyNE-based M1 model.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim




# -----------------------------------------------------------
# Main code
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  # create network object and set cfg and net params

sim.pc.timeout(300)                          # set nrn_timeout threshold to X sec (max time allowed without increasing simulation time, t; 0 = turn off)
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
#sim.net.connectCells()            			# create connections between cells based on params
#sim.net.addStims() 							# add network stimulation
#sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)

# Simulation option 1: standard
#sim.runSim()                              # run parallel Neuron simulation (calling func to modify mechs)

#print(cfg.modifyMechs)
# Simulation option 2: interval function to modify mechanism params
#sim.runSimWithIntervalFunc(1000.0, modifyMechsFunc)       # run parallel Neuron simulation (calling func to modify mechs)

# Gather/save data option 1: standard
#sim.gatherData()

# Gather/save data option 2: distributed saving across nodes 
#sim.saveDataInNodes()
#sim.gatherDataFromFiles()

#sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
#sim.analysis.plotData()         			# plot spike raster etc

