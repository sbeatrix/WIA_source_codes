## This script generates the stochastic trajectories necessary to create
## figures for the stochastic version of the model simulating cells where the drug was washed out after x minutes.
## WT, APC-A, GM2

# Import necessary libraries
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for matplotlib
import os  # For folder creation in code
import pickle  # For saving large result files
import stochpy
smod = stochpy.SSA()  # Initialize stochastic simulation algorithm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm  # Colormap for plotting

# Set figure resolution
plt.rcParams['figure.dpi'] = 150

# Import stochastic model
smod.Model('model_washout.psc', dir='../models_psc/')

### Simulation settings
traj_nr = 100  # Number of trajectories to simulate
time = 400  # Simulation end time

# Simulating wild-type (WT) washout
res_wo_wt = []  # List to store simulation results for wild-type
for n in np.arange(1, traj_nr+1):
    smod.DoStochSim(method="direct", end=time, mode="time", species_selection=['ac', 'nuk'])
    res_wo_wt.append(smod.data_stochsim.getSpecies())

# Save WT washout simulation results
f=open('../trajectories/res_wt_wo.pickle', 'w')
pickle.dump(res_wo_wt, f)
f.close()

#####################################################################################################################

# Simulating APC-A mutant washout
new_kassac = 0.07  
smod.ChangeParameter("kassac", new_kassac)

# Set initial conditions for APC-A mutant (determined using XPP)
initial_conditions_apca = {
    "a": 56, "ac": 14, "acmc": 29, "mc": 7, "mstar": 137,
    "m": 1, "c": 14, "x": 28
}
for var, val in initial_conditions_apca.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

res_wo_apca = []  # List to store simulation results for APC-A mutant
for n in np.arange(1, traj_nr+1):
    smod.DoStochSim(method="direct", end=time, mode="time", species_selection=['ac', 'nuk'])
    res_wo_apca.append(smod.data_stochsim.getSpecies())

# Save APC-A mutant washout simulation results
f=open('../trajectories/res_apca_wo.pickle', 'w')
pickle.dump(res_wo_apca, f)
f.close()

# Restore original kassac value
orig_kassac = 0.1073
smod.ChangeParameter("kassac", orig_kassac)

######################################################################################################################

# Simulating GM2 mutant washout
# Set initial conditions for GM2 mutant (determined using XPP)
initial_conditions_gm2 = {
    "a": 52, "ac": 18, "acmc": 29, "mc": 6, "mstar": 183,
    "m": 2, "c": 10, "x": 23
}
for var, val in initial_conditions_gm2.items():
    smod.ChangeInitialSpeciesCopyNumber(var, val)

res_wo_gm2 = []  # List to store simulation results for GM2 mutant
for n in np.arange(1, traj_nr+1):
    smod.DoStochSim(method="direct", end=time, mode="time", species_selection=['ac', 'nuk'])
    res_wo_gm2.append(smod.data_stochsim.getSpecies())

# Save GM2 mutant washout simulation results
f=open('../trajectories/res_gm2_wo.pickle', 'w')
pickle.dump(res_wo_gm2, f)
f.close()