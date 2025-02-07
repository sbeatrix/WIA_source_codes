### Parameter fitting using the deterministic XML model
# This script is used for parameter fitting. We defined constraints and reasonable intervals
# based on experiments and literature.  It also includes optimization
# using the Levenberg-Marquardt algorithm and parameter ensemble analysis.

# Import necessary libraries
import matplotlib
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import SCfunctions as sc
import libsbml as lib
import numpy as np

# Set figure resolution
plt.rcParams['figure.dpi'] = 100

# Load the deterministic model from an SBML file
net_wt = IO.from_SBML_file("../models_det/model_det_ext.xml", "net_wt")

# Compile the model
try:
    net_wt.compile()
except:
    net_wt.compile()

## Define initial parameter set (has to be bistable)
parameter_values = {
    "ksyn_c20": 37,
    "ksyn_cy": 320,
    "kdegbg": 0.023,
    "kinact": 0.45,
    "J": 0.062,
    "kact": 0.0014,
    "k_nuk": 78,
    "J1": 0.50,
    "kassmc": 0.022,
    "kdissmc": 0.63,
    "kassac": 0.023,
    "kassacmc": 0.045,
    "kdissac": 0.66,
    "kdissacmc": 0.43,
    "kdeg": 1.1,
    "kdeg_cy": 0.024,
    "kdegbg_cy": 0.068,
    "nuk": 32
}

# Set initial conditions for biochemical species
initial_conditions = {
    "C": 25, "Mstar": 106, "MC": 38, "AC": 30, "ACMC": 30,
    "M": 1, "A": 80, "CycB": 430
}

# Total amounts of molecular components
totals = {"Ctot":150, "Mtot":175, "Atot":150}

# Apply parameter and initial conditions to the network
net_wt.set_var_ics({**parameter_values, **initial_conditions, **totals})

# Fix certain parameters as non-optimizable
for param in ["Mtot", "Atot", "nuk"]:
    net_wt.set_var_optimizable(param, False)

# Define experimental data for validation
exp_data = {
    "net_wt": {
        "MC": {10: (36, 1), 1000: (36, 1)},
        "ACMC": {10: (30, 0.1), 1000: (30, 0.1)},
        "A": {10: (60, 10), 1000: (60, 10)},
        "CycB": {10: (400, 10), 1000: (400, 10)},
        "Ctot": {10: (150, 5), 1000: (150, 5)}
    }
}
exp_concentrations = Experiment(name="concentrations", data=exp_data)
exp_concentrations.set_fixed_sf({"MC": 1, "ACMC": 1, "A": 1, "CycB": 1, "Ctot": 1})

## Define initial conditions for SAC on and SAC off states
var_ics_list = [
    {"C": 0, "MC": 0, "AC": 0, "ACMC": 0, "CycB": 0, "Mstar": 180, "nuk": 32},
    {"C": 0, "MC": 0, "AC": 150, "ACMC": 0, "CycB": 0, "Mstar": 0, "nuk": 32},
    {"C": 0, "MC": 0, "AC": 0, "ACMC": 0, "CycB": 0, "Mstar": 180, "nuk": 1},
    {"C": 0, "MC": 0, "AC": 150, "ACMC": 0, "CycB": 0, "Mstar": 0, "nuk": 1}
]

# Copy network and set initial conditions for different states
nets = [net_wt.copy(new_id=f"net_{i+1}") for i in range(4)]
trajectories = []
for i, net in enumerate(nets):
    net.set_var_ics(var_ics_list[i])
    trajectories.append(net.integrate([0, 2000]))

# Define experimental data to enforce stable steady states
exp_data_bistability = {
    f"net_{i+1}": {
        "AC": {
            500: (10 if i % 2 == 0 else 100, 10 if i % 2 == 0 else 20),
            1500: (10 if i % 2 == 0 else 100, 10 if i % 2 == 0 else 20),
        }
    }
    for i in range(4)
}
exp_bistability = Experiment(name='bistability', data=exp_data_bistability)
exp_bistability.set_fixed_sf({"AC": 1.0})

# Define model using experiments and networks
m = Model([exp_concentrations, exp_bistability], [net_wt] + nets)
params = m.get_params()

# Add residuals for prior information in log space
for var, val in params.items():
    res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(10)))
    m.AddResidual(res)

# Perform parameter optimization
params_opt = Optimization.fmin_lm_log_params(m, params, maxiter=1000, disp=True)

# Run ensemble analysis
Network.full_speed()
j = m.jacobian_log_params_sens(np.log(params_opt))
jtj = np.dot(np.transpose(j), j)
ens, gs, r = Ensembles.ensemble_log_params(m, params_opt, jtj, steps=100000, temperature=1)
ens_pruned = ens[::1000]

# Plot parameter distributions
n_rows, n_cols = 2, 7
fig, axes = plt.subplots(n_rows, n_cols, figsize=(11,3), sharey=True)
logbins = np.logspace(np.log10(1e-2), np.log10(1e2), 30)
for i, ax in enumerate(axes.flat):
    if i == len(params_opt):
        break
    var = list(params_opt.keys())[i]
    val = np.array([e[i] for e in ens_pruned]) / params.getByKey(var)
    ax.set_xscale('log')
    ax.hist(val, color="white", edgecolor="black", bins=logbins)
    ax.set_title(var)
plt.tight_layout()
fig.savefig("../ensemble.svg")
