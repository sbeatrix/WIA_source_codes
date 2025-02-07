# WIA_source_codes

## model_files/
This folder contains all the files required for running the provided scripts and bifurcation analysis in XPP-AUTO.

### deterministic/model_det_ext.xml
is an xml file describing the deterministic model. It is used by parameter_fitting.py script (see in ../scripts). It is converted to SBML file.

### deterministic/model_det.ode
is the .ode file required for steady state and bifurcation analysis in XPP of the complex model.

### deterministic/simplemodel1.ode
is the .ode file required for steady state and bifurcation analysis in XPP of the first simple model.

### deterministic/simplemodel2.ode
is the .ode file required for steady state and bifurcation analysis in XPP of the second simple model.

### stochastic/model_arrest.psc 
is a .psc file required by the stochastic package. It is used for simulating cells under constant nocodazole arrest in the complex model.

### stochastic/model_washout.psc
is a .psc file required by the stochastic package. It is used for simulating cells in condition where the drug is washed out in the complex model.

### stochastic/model_arrest_s1.psc 
is a .psc file required by the stochastic package. It is used for simulating cells under constant nocodazole arrest in the first simple model.

### stochastic/model_washout_s1.psc
is a .psc file required by the stochastic package. It is used for simulating cells in condition where the drug is washed out in the first simple model.

### stochastic/model_arrest_s2.psc 
is a .psc file required by the stochastic package. It is used for simulating cells under constant nocodazole arrest in the second simple model.

### stochastic/model_washout_s2.psc
is a .psc file required by the stochastic package. It is used for simulating cells in condition where the drug is washed out in the second simple model.


## scripts/
This folder contains the scripts used to fit the model to experimental data, for studying the behavior of the model and to perform stochastic simulations to reproduce experiments.

### parameter_fitting.py
This script uses the deterministic model, introduces the constraints defined, based on measurements and experimental data. First it looks for the local minimum, then starting from the optimized parameter set it builds a parameter ensemble which consists of 100 individual parameter sets which fulfill the given constraints.

### Stochastic simulations:

#### Input for simulations:
For each stochastic simulation a .psc model file is needed, describing the reactions of the system. There is a separate model used for constant drug arrest and washout simulations (see ../model_files/stochastic)

#### Simulate trajectories
Scripts **scripts/maketraj_arrest.py** and **scripst/maketraj_washout.py** were used to run the stochastic simulations 
for wildtype and both mutants in case of constant drug arrest and drug washout. These scripts save out trajectories to .pickle files, used later for processing data.

#### Process data for figures
**scripts/makefigures.py** processes the saved trajectories and generates data for ‘simulation’ figures: Fig 2D-E, 3A-C, 4B-D, 5B-D, S3B-C, S4-C, S5-A

#### Simple models
For simple models the same scripts can be used to generate subfigures of FigS6, changing the input files accordingly to: stochastic/model_arrest_s1.psc, stochastic/model_arrest_s2.psc, stochastic/model_washout_s1.psc, stochastic/model_washout_s2.psc.

