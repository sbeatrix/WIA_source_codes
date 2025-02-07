### Import and Configure Libraries
import matplotlib
matplotlib.use('Agg')  # Set backend to 'Agg' for non-interactive environments
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import csv

# Configure Matplotlib settings for consistent figure rendering
plt.rcParams.update({
    'figure.dpi': 100,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

## Functions

def switch_nuk(traj, species_col, threshold):
    """ 
    Determine the timepoint of adaptation and the number of unattached kinetochores at adaptation. 
    
    Parameters:
        traj (array): Stochastic simulation result trajectories
        species_col (int): Column index of the species from which the threshold is calculated
        threshold (int): Defined threshold for transition from SAC-on to SAC-off

    Returns:
        list: [Timepoint of adaptation, Number of unattached kinetochores at adaptation]
        If adaptation does not occur, returns [infinity, maximum unattached kinetochores].
    """
    t_switch_up = np.argmax(traj[:, species_col] > threshold)
    return [traj[t_switch_up, 0], traj[t_switch_up, 2]] if t_switch_up != 0 else [np.inf, 10]
   
def pad_or_truncate(trajs, target_length):
    """ 
    Pad or truncate trajectories to a uniform length.
    
    Parameters:
        trajs (list): List of trajectories
        target_length (int): Minimum required length for all trajectories

    Returns:
        np.array: Array of padded or truncated trajectories
    """
    return np.array([
        np.vstack([traj, np.tile(traj[-1, :], (target_length - len(traj), 1))])
        if len(traj) < target_length else traj[:target_length]
        for traj in trajs
    ])

## Load Simulated Trajectories

def load_pickle(filepath):
    """ Load a pickle file and return the data """
    with open(filepath, 'rb') as file:
        return pickle.load(file)

res_wt_arrest = load_pickle('../trajectories/res_wt_arrest.pickle')
res_apca_arrest = load_pickle('../trajectories/res_apca_arrest.pickle')
res_gm2_arrest = load_pickle('../trajectories/res_gm2_arrest.pickle')
res_wt_wo = load_pickle('../trajectories/res_wt_wo.pickle')
res_apca_wo = load_pickle('../trajectories/res_apca_wo.pickle')
res_gm2_wo = load_pickle('../trajectories/res_gm2_wo.pickle')

## Compute Adaptation Times for Drug Arrest Conditions
threshold = 80  # Define threshold based on steady states
switches_wildtype, switches_apca, switches_gm2 = [], [], []

# Compute adaptation times for different mutant types
for r in res_wt_arrest:
    switches_wildtype.append(switch_nuk(r, 1, threshold)[0])
for r in res_apca_arrest:
    switches_apca.append(switch_nuk(r, 1, threshold)[0])
for r in res_gm2_arrest:
    switches_gm2.append(switch_nuk(r, 1, threshold)[0])

# Sort adaptation times
switches_wildtype.sort()
switches_apca.sort()
switches_gm2.sort()

## Save Adaptation Times for ECDF Plots

def save_to_csv(filename, data):
    """ Save adaptation times to CSV file """
    with open(filename, mode='wb') as file:
        writer = csv.writer(file)
        writer.writerow(data)

save_to_csv('../simulation_results/fullmodel_wt_arrest_ecdf_data.csv', switches_wildtype)
save_to_csv('../simulation_results/fullmodel_apca_arrest_ecdf_data.csv', switches_apca)
save_to_csv('../simulation_results/fullmodel_gm2_arrest_ecdf_data.csv', switches_gm2)

## Save Washout Simulation Results

def process_washout_data(res_data, filename):
    """ Compute and save washout simulation results """
    switches = [switch_nuk(i, 1, threshold) for i in res_data]
    with open(filename, mode='wb') as file:
        writer = csv.writer(file)
        writer.writerow(['Time', 'nUK'])
        writer.writerows(switches)

process_washout_data(res_wt, '../simulation_results/washout_wt_nuk_at_adapt.csv')
process_washout_data(res_apca, '../simulation_results/washout_apca_nuk_at_adapt.csv')
process_washout_data(res_gm2, '../simulation_results/washout_gm2_nuk_at_adapt.csv')

## Compare Kinetochore Attachment Dynamics: Simulation vs. Experiment

# Determine the minimum trajectory length
min_length = min(traj.shape[0] for traj in res_wt)
res_wt_padded = pad_or_truncate(res_wt, min_length)
n_time_points = res_wt_padded[0].shape[0]

# Initialize arrays for mean and standard deviation
means_nuk = np.zeros(n_time_points)
stds_nuk = np.zeros(n_time_points)

# Compute mean and standard deviation for each time point
for t in range(n_time_points):
    nuk_values = np.array([traj[t, 2] for traj in res_wt_padded])
    means_nuk[t] = np.mean(nuk_values)
    stds_nuk[t] = np.std(nuk_values)

# Load experimental data
file_path = '../exp_data/WIA_Overall_Mad2_delocalization.txt'
df = pd.read_csv(file_path, delimiter='\t')
df = df.iloc[17:].reset_index(drop=True)

# Normalize experimental data
last_ten_avg = df["fracSACON"].iloc[-10:].mean()
adjusted_fracSACON = df["fracSACON"] - last_ten_avg
rescaled_fracSACON = (adjusted_fracSACON - adjusted_fracSACON.min()) / (adjusted_fracSACON.max() - adjusted_fracSACON.min())

# Plot simulation vs. experimental results
fig, ax = plt.subplots()
plt.plot(df["Time"], rescaled_fracSACON, marker='o', color='grey', alpha=0.8, label='Experimental')
plt.plot(time_points-wo, means_nuk/10, label='Simulation', color='k')
plt.fill_between(time_points-wo, (means_nuk - stds_nuk)/10, (means_nuk + stds_nuk)/10, color='grey', alpha=0.2)
plt.xlabel("Time after washout (min)", fontsize=12)
plt.ylabel("fracSACON", fontsize=12)
plt.xlim([-3, 302])
plt.legend()
fig.savefig('../simulation_results/figures/FigS3B.pdf')
plt.show()
