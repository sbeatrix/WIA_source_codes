### Import and Configure Libraries
import matplotlib
matplotlib.use('Agg') 
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import csv
# Configure Matplotlib
plt.rcParams.update({
    'figure.dpi': 100,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})


## Functions

# switch_nuk: function to define the number of unattached kineotchores and the time point the SAC on - SAC off switch happens. 
# We follow AC transition from the lower (SAC on) to upper steady state (SAC off)
# input:
## traj - stochastic simulation result trajectories
## species_col - the column index of the species from species selection, from the steady states of which we calculate the threshold
# output[0] - timepoint of adaptation
# output[1] - number of unattached kinetochores at the time of adaptation
# In case a trajectory never adapts in the simulated time, the function return infinity and maximum nr of unattached kinetochores

def switch_nuk(traj, species_col, threshold):
    """ Determine the timepoint of adaptation and number of unattached kinetochores at adaptation. """
    t_switch_up = np.argmax(traj[:, species_col] > threshold)
    return [traj[t_switch_up, 0], traj[t_switch_up, 2]] if t_switch_up != 0 else [np.inf, 10]
   
# pad_or_truncate: function needed for plotting the nuk attachment dynamics mean and deviation of the trajectories    
# pad or truncate all trajectories to the same length
# input: [0] trajs - trajectories, here nuk trajectories; [1] target_length - give the minimum length 
# output - padded trajectories
def pad_or_truncate(trajs, target_length):
    """ Pad or truncate trajectories to a uniform length. """
    return np.array([
        np.vstack([traj, np.tile(traj[-1, :], (target_length - len(traj), 1))])
        if len(traj) < target_length else traj[:target_length]
        for traj in trajs
    ])

# open the file containing the simulated trajectories for drug arrest, wild-type
file = open('../trajectories/res_wt_arrest.pickle', 'rb')
res_wt_arrest= pickle.load(file)
file.close()

# open the file containing the simulated trajectories for drug arrest, APC-A mutant
file_apca = open('../trajectories/res_apca_arrest.pickle','rb')
res_apca_arrest= pickle.load(file_apca)
file_apca.close()


# open the file containing the simulated trajectories for drug arrest, GM2 mutant
file_gm2 = open('../trajectories/res_gm2_arrest.pickle','rb')
res_gm2_arrest= pickle.load(file_gm2)
file_gm2.close()

# open the file containing the simulated trajectories for drug washout, wild-type
file = open('../trajectories/res_wt_wo.pickle', 'rb')
res_wt_wo= pickle.load(file)
file.close()

# open the file containing the simulated trajectories for drug washout, APC-A
file = open('../trajectories/res_apca_wo.pickle', 'rb')
res_apca_wo= pickle.load(file)
file.close()

# open the file containing the simulated trajectories for drug washout,GM2 mutant
file = open('../trajectories/res_gm2_wo.pickle', 'rb')
res_gm2_wo= pickle.load(file)
file.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
## Save data necessary for the cumulative distribution for adaptation times, drug arrest (constant nuk), wild-type
threshold = 80 # threshold defined based on the steady states
switches_wildtype = list()
switches_apca = list()
switches_gm2 = list()

# iterate through wild-type, apc-a, gm2 arrest trajectories and save the timepoints of adaptations to a list
for r in res_wt_arrest:
    switches_wildtype.append(switch_nuk(r,1,threshold)[0])
    
for r in res_apca_arrest:
    switches_apca.append(switch_nuk(r,1,threshold)[0])  

for r in res_gm2_arrest:
    switches_gm2.append(switch_nuk(r,1,threshold)[0])
    
switches_wildtype.sort()  
switches_apca.sort()  
switches_gm2.sort()  

# save the data for creating ECDF (used for Fig 2E)
filename='../simulation_results/fullmodel_wt_arrest_ecdf_data.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(switches_wildtype)
    
# save the data for creating ECDF (used for Fig S4-C)
filename='../simulation_results/fullmodel_apca_arrest_ecdf_data.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(switches_apca)

# save the data for creating ECDF (used for Fig S5-A)
filename='../simulation_results/fullmodel_wt_arrest_ecdf_data.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(switches_gm2)
  
#---------------------------------------------------------------------------------------------------------------------------------------
## save out for washout simulations time of adaptation and nuk at adaptation, for wt, apc-a and gm2
## creates input data to generate Figures: 3A, 3C; 4B, 4C, 4D; 5B, 5C, 5D

switches_wt=list()
for i in res_wt:
    switches_wt.append(switch_nuk(i,1,threshold))
    # Write data to CSV file
filename='../simulation_results/washout_wt_nuk_at_adapt.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(['Time', 'nUK'])
    writer.writerows(switches_wt)
    
switches_apca=list()
for i in res_apca:
    switches_apca.append(switch_nuk(i,1,threshold))
    # Write data to CSV file
filename='../simulation_results/washout_apca_nuk_at_adapt.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(['Time', 'nUK'])
    writer.writerows(switches_apca)

switches_gm2=list()
for i in res_gm2:
    switches_gm2.append(switch_nuk(i,1,threshold))
    # Write data to CSV file
filename='../simulation_results/washout_gm2_nuk_at_adapt.csv'
with open(filename, mode='wb') as file:
    writer = csv.writer(file)
    writer.writerow(['Time', 'nUK'])
    writer.writerows(switches_gm2)
    

#---------------------------------------------------------------------------------------------------------------------------------------
### Compare kinetochore attachment dynamics for simulation and experiments (Fig S3B)

min_length = min(traj.shape[0] for traj in res_wt)

res_wt_padded = pad_or_truncate(res_wt, min_length)

# Check the shape of the padded trajectories
n_time_points = res_wt_padded[0].shape[0]

# Initialize arrays which contain mean and standard deviation values
means_nuk = np.zeros(n_time_points)
stds_nuk = np.zeros(n_time_points)

# Calculate the mean and standard deviation for the trajectories
for t in range(n_time_points):
    nuk_values = np.array([traj[t, 2] for traj in res_wt_padded])
    means_nuk[t] = np.mean(nuk_values)
    stds_nuk[t] = np.std(nuk_values)

# Plot the mean and the extent of fluctuations
time_points = res_wt_padded[0][:, 0]

## read the experimental data
file_path = '../exp_data/WIA_Overall_Mad2_delocalization.txt'  
df = pd.read_csv(file_path, delimiter='\t')

# Ignore values before the washout
df = df.iloc[17:].reset_index(drop=True)

# Calculate the average of the last ten points
last_ten_avg = df["fracSACON"].iloc[-10:].mean()

# Subtract the average of the last ten points from all points
adjusted_fracSACON = df["fracSACON"] - last_ten_avg

# Rescale between 0 and 1
min_val = adjusted_fracSACON.min()
max_val = adjusted_fracSACON.max()
rescaled_fracSACON = (adjusted_fracSACON - min_val) / (max_val - min_val)
fig, ax=plt.subplots()

# Plot the rescaled data
plt.plot(df["Time"], rescaled_fracSACON, marker='o',color='grey', alpha=0.8, label='Experimental')


# Plot nuk trajectory mean and deviation
plt.plot(time_points-wo, means_nuk/10, label='Simulation', color='k')
plt.fill_between(time_points-wo, (means_nuk - stds_nuk)/10, (means_nuk + stds_nuk)/10, color='grey', alpha=0.2)
plt.xlabel("Time after washout (min)", fontsize=12)
plt.ylabel("fracSACON", fontsize=12)
plt.xlim([-3,302])
plt.legend()
fig.savefig('../simulation_results/figures/FigS3B.pdf')
plt.show()
