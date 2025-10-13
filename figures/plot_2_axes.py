import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Serif"
plt.rcParams["font.size"] = 17.0
plt.rcParams["axes.labelsize"] = 17.0
plt.rcParams['lines.linewidth'] = 1.0 # 1.0
plt.rcParams["xtick.minor.visible"] = True 
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"]= True 
plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] =True
plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 0.75
plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 0.75
plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 5.0
plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 2.5
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams['figure.dpi'] = 300.0
## change legend parameters
plt.rcParams["legend.fontsize"] = 17.0
plt.rcParams["legend.frameon"] = True
subdict = {"figsize" : (3.25,3.5),"constrained_layout" : True,"sharex" : True}

# === Settings ===
file_location = "figures/arrow/1_arrow_output.txt"  # update with your actual filename
x_is_log_scale = False
y_is_log_scale = False


# Axis control
x_limit = (0.0, 1.4)
x_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
x_tick_labels = ["0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4"]

y_limit_left = (175.0, 215.0)  # velocity axis
y_ticks_left = [175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0]
y_tick_labels_left = ["175", "180", "185", "190", "195", "200", "205", "210", "215"]

y_limit_right = (-0.015, 0.025)  # angle of attack axis
y_ticks_right = [-0.015,-0.01,-0.005,0.000,0.005,0.010,0.015,0.020,0.025]
# y_ticks_right = [-1.0,1.0]
y_tick_labels_right = ["-0.015","-0.01","-0.005","0.000","0.005","0.010","0.015","0.020","0.025"]
# y_tick_labels_right = ["-1.0","1.0"]

# === Data Reader ===
def read_csv_data(file_path):
    """
    Reads a CSV file, ignoring commented lines (#).
    Returns numpy array of shape (n, m).
    """
    return np.loadtxt(file_path, skiprows = 1)

# === Load Data ===
data = read_csv_data(file_location)

# Columns: adjust if needed
time = data[:, 0]       # column 0
velocity_x = data[:, 1]   # column 1
velocity_y = data[:, 2] 
velocity_z = data[:, 3] 
velocity = np.sqrt(velocity_x**2 + velocity_y**2 + velocity_z**2)
x_comp = data[:, 1]     # column 1 (used for atan2 denominator)
y_comp = data[:, 3]     # column 3 (used for atan2 numerator)

# Compute angle of attack in degrees
alpha = np.arctan2(y_comp, x_comp)

# === Plotting ===
fig, ax1 = plt.subplots()

# Left y-axis: velocity
ax1.plot(time, velocity, 'k-', label="V [ft/s]")
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Velocity [ft/s]", color="k")
ax1.tick_params(axis="y", labelcolor="k")

# Apply x-axis formatting
ax1.set_xlim(x_limit)
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_tick_labels)

# Apply left y-axis formatting
ax1.set_ylim(y_limit_left)
ax1.set_yticks(y_ticks_left)
ax1.set_yticklabels(y_tick_labels_left)

# Right y-axis: angle of attack
ax2 = ax1.twinx()
ax2.plot(time, alpha, 'k--', label="$\\alpha$ [rad]")
ax2.set_ylabel("Angle of Attack [rad]", color="k")
ax2.tick_params(axis="y", labelcolor="k")

# Apply right y-axis formatting
ax2.set_ylim(y_limit_right)
ax2.set_yticks(y_ticks_right)
ax2.set_yticklabels(y_tick_labels_right)

# Optional scaling
if x_is_log_scale:
    ax1.set_xscale("log")
if y_is_log_scale:
    ax1.set_yscale("log")
    ax2.set_yscale("log")

# Merge legends from both axes
lines, labels_left = ax1.get_legend_handles_labels()
lines2, labels_right = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels_left + labels_right, loc="upper right")

plt.title("Velocity and Angle of Attack vs Time")
plt.tight_layout()
ax1.grid(True)
ax2.grid(True)
plt.savefig("figures/arrow/alpha_plot.svg")