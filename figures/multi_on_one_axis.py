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
file_location = "figures/arrow/1_arrow_output.txt"

# Axis control
x_limit = (0.0, 1.4)
x_ticks = np.arange(0.0, 1.41, 0.2)
x_tick_labels = [f"{val:.1f}" for val in x_ticks]

y_limit = (-1.5, 1.5)
y_ticks = np.arange(-1.5, 1.51, 0.5)
y_tick_labels = [f"{val:.1f}" for val in y_ticks]

# === Data Reader ===
def read_txt_data(file_path):
    """
    Reads a whitespace-delimited text file with one header row.
    """
    return np.loadtxt(file_path, comments="#", skiprows=1)

# === Load Data ===
data = read_txt_data(file_location)
x_comp = data[:, 1]     # column 1 (used for atan2 denominator)
v_comp = data[:, 2]
y_comp = data[:, 3]     # column 3 (used for atan2 numerator)
velocity = np.sqrt(x_comp**2 + v_comp**2 + y_comp**2)
# Compute angle of attack in degrees
alpha = np.rad2deg(np.arctan2(y_comp, x_comp))
beta = np.rad2deg(np.arcsin(v_comp/velocity))

# Columns (adjust as needed)
x_col = data[:, 0]   # shared x-axis
y1_col = alpha  # first y column
y2_col = beta  # second y column

# === Plotting ===
fig, ax = plt.subplots()

# Plot both y-series on the same axis
ax.plot(x_col, y1_col, 'k-', label="$\\alpha [deg]$")
ax.plot(x_col, y2_col, 'k--', label="$\\beta [deg]$")

# Axis labels
ax.set_xlabel("Time (s)")
ax.set_ylabel("Angle [deg]")

# Apply formatting
ax.set_xlim(x_limit)
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_tick_labels)

ax.set_ylim(y_limit)
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_tick_labels)

# Legend
ax.legend(loc="upper right")

plt.title("$\\alpha$ and $\\beta$ vs Time")
plt.savefig("figures/arrow/alpha_vs_beta_plot.svg")
plt.tight_layout()