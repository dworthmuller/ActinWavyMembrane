import pandas as pd
import numpy as np
from scipy import interpolate          # USED
from mpl_toolkits.axes_grid1 import make_axes_locatable  # UNUSED
import matplotlib.pyplot as plt        
from matplotlib.lines import Line2D    # USED
import proplot as pplt                 # USED
from scipy import integrate            # USED (cumtrapz)
import glob                            # USED

## SET SOME GENERAL PARAMS AND LOAD COLORS YOU WANT TO USE ###
cc = ['#F3DECA','#FA9483','#2D4057','#4097aa','#453f3d']  # UNUSED
cc_special = ['#5eac82']                                  # UNUSED

# Plot style dictionaries (currently not passed to pplt.rc or context → UNUSED)
params = {
    'subplots.share': False,
    'subplots.span': False,
    'metawidth': 1.5,
    'gridwidth': 0,
    'legend.fontsize': 8,
    'font.sans-serif': 'Source Sans Pro'
}  # UNUSED; also overwritten by 'params = df.iloc[0].to_dict()' below

params_dark = {                     # UNUSED
    'subplots.share': False,
    'subplots.span': False,
    'metawidth': 1.5,
    'gridwidth': 0,
    'legend.fontsize': 8,
    'font.sans-serif': 'Source Sans Pro',
    "lines.color": "white",
    "patch.edgecolor": "white",
    "text.color": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "white",
    "axes.labelcolor": "white",
    "xtick.color": "white",
    "xtick.labelcolor": "white",
    "ytick.labelcolor": "white",
    "ytick.color": "white",
    "grid.color": "white",
    "axes.facecolor": "black",
    "figure.facecolor": "black",
    "figure.edgecolor": "black",
    "savefig.facecolor": "black",
    "savefig.edgecolor": "black",
    "boxplot.boxprops.color": "white",
    "boxplot.capprops.color": "white",
    "boxplot.flierprops.color": "white",
    "boxplot.flierprops.markeredgecolor": 'k'"white",  # This looks like a typo: 'k'"white"
    "boxplot.whiskerprops.color": "white"
}

# Colormaps defined but not used in this script
vik = pplt.Colormap('vik')       # UNUSED
bilbao = pplt.Colormap('bilbao') # UNUSED

# -------------------------------------------------------------------
# Paths and parameters
# -------------------------------------------------------------------
__DATAPATH__ = 'Data/'
__SAVEPATH__ = './'

# Read CSV containing simulation parameters and previously saved forces
df = pd.read_csv('inputParams_add.csv')
params = df.iloc[0].to_dict()   # NOTE: overwrites the 'params' plotting dict defined above

# Count how many surface-stress CSV files we have
no_csv_files = len(glob.glob('./Data/surface_stresses*.{}'.format('csv')))

# Time (or iteration) index array based on the number of CSVs
iteration_steps = np.arange(0, no_csv_files, 1)

# Containers: stress fields as functions of x and time
sig_yy_x_vs_t = []
sig_nn_x_vs_t = []
f_y_x_vs_t = []
x_vs_t = []
shape_deriv_x_vs_t = []

# -------------------------------------------------------------------
# Load all surface stress profiles into memory
# -------------------------------------------------------------------
for step in iteration_steps:
    data = pd.read_csv(__DATAPATH__ + "surface_stresses%04d.csv" % (step))
    x = data["x"].to_numpy()
    sig_nn = data["sig_nn"].to_numpy()
    sig_yy = data["sig_yy"].to_numpy()
    f_y = data["f_y"].to_numpy()
    shape_deriv = data["shape_deriv"].to_numpy()

    sig_nn_x_vs_t.append(sig_nn)
    sig_yy_x_vs_t.append(sig_yy)
    f_y_x_vs_t.append(f_y)
    shape_deriv_x_vs_t.append(shape_deriv)
    x_vs_t.append(x)

# -------------------------------------------------------------------
# STANDARD LINE PLOT: σ_nn for all time steps + final T_y and σ_yy
# -------------------------------------------------------------------
fig = pplt.figure(share=False, aspect=1.61)
ax1 = fig.subplot(111)

# Color cycles: one color per iteration step
cycle_reds = pplt.Cycle('reds', no_csv_files, left=0.1)
cycle_blues = pplt.Cycle('blues', no_csv_files, left=0.1)
cycle_greens = pplt.Cycle('greens', no_csv_files, left=0.1)
cc_custom_reds = cycle_reds.by_key()['color']
cc_custom_blues = cycle_blues.by_key()['color']
cc_custom_greens = cycle_greens.by_key()['color']

# Common fine x-grid for interpolation (based on first profile)
xnew = np.arange(x_vs_t[0][0], x_vs_t[0][-1] + 0.01, 0.01)

# Plot σ_nn for all time steps
for step in iteration_steps:
    x = x_vs_t[step]
    f_sig_nn = interpolate.interp1d(x, sig_nn_x_vs_t[step], kind='linear')
    ax1.plot(xnew, f_sig_nn(xnew), color=cc_custom_blues[step], label=r'$\sigma_{nn}$')

# For final time step: plot σ_yy and T_y
# NOTE: 'x' here is still defined from the last loop iteration (i.e. x_vs_t[-1])
x = x_vs_t[-1]  # optional: make this explicit
f_sig_yy = interpolate.interp1d(x, sig_yy_x_vs_t[-1], kind='linear')
f_f_y = interpolate.interp1d(x, f_y_x_vs_t[-1], kind='linear')
ax1.plot(xnew, f_f_y(xnew), color=cc_custom_greens[-1], label=r'$T_{y}$')
ax1.plot(xnew, f_sig_yy(xnew), color=cc_custom_reds[-1], label=r'$\sigma_{yy}$')

# Custom legend using representative colors
legend_elements = [
    Line2D([0], [0], color=cc_custom_blues[-1], lw=2, label=r'$\sigma_{nn}$'),
    Line2D([0], [0], color=cc_custom_greens[-1], lw=2, label=r'$T_{y}$'),
    Line2D([0], [0], color=cc_custom_reds[-1], lw=2, label=r'$\sigma_{yy}$')
]
ax1.legend(handles=legend_elements, loc='t', ncols=3, frame=False)

# Clean up spines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Axis labels and limits
ax1.format(
    xlabel=r'x',
    ylabel=r'stress (force)',
    xlim=(x_vs_t[0][0], x_vs_t[0][-1])
)

fig.save(__SAVEPATH__ + 'surface_stresses.pdf')
fig.save(__SAVEPATH__ + 'surface_stresses.png')

# -------------------------------------------------------------------
# INTEGRAL OF Y-TRACTION: ∫ T_y sqrt(1+f'^2) dx for each time step
# -------------------------------------------------------------------
fig = pplt.figure(share=False, aspect=1.61)
ax1 = fig.subplot(111)

cycle_grays = pplt.Cycle('grays', no_csv_files, left=0.1)
cc_custom_grays = cycle_grays.by_key()['color']

# xnew defined again but not used in this block
xnew = np.arange(x_vs_t[0][0], x_vs_t[0][-1] + 0.01, 0.01)  # UNUSED below

for step in iteration_steps:
    x = x_vs_t[step]
    x = x_vs_t[-1]   # NOTE: this overwrites the previous line; the first assignment is redundant

    # Local interpolant for T_y at this step (not used in final code)
    f_f_y = interpolate.interp1d(x, f_y_x_vs_t[step], kind='cubic')  # UNUSED; integration uses raw array

    # Arc-length–weighted integral of T_y along the interface:
    # ∫ T_y sqrt(1+f'^2) dx
    y_int = integrate.cumtrapz(
        f_y_x_vs_t[step] * np.sqrt(1 + shape_deriv_x_vs_t[step]**2),
        x,
        initial=0
    )
    # Alternative with interpolant on fine grid (commented out)
    # y_int = integrate.cumtrapz(f_f_y(xnew), xnew, initial=0)

    ax1.plot(x, y_int, color=cc_custom_grays[step], label=r'$\int\;T_y dx$')

# Reference lines
ax1.axhline(0, linestyle='--', color='black', lw=1)
ax1.axvline(2, linestyle='--', color='red', lw=1)

legend_elements = [
    Line2D([0], [0], color=cc_custom_grays[-1], lw=2, label=r'$\int\;T_y dx$')
]
ax1.legend(handles=legend_elements, loc='t', ncols=3, frame=False)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.format(
    xlabel=r'x',
    ylabel=r'stress (force)',
    xlim=(x_vs_t[0][0], x_vs_t[0][-1])
)

fig.save(__SAVEPATH__ + 'integral_T_y.pdf')
fig.save(__SAVEPATH__ + 'integral_T_y.png')

# -------------------------------------------------------------------
# Save steady-state (final time) stress values at specific positions
# -------------------------------------------------------------------
# NOTE: here 'x' and 'step' refer to their last values from the previous loop:
#       x = x_vs_t[-1], step = last element of iteration_steps.

# Use last profile as steady state
x = x_vs_t[-1]
f_sig_nn = interpolate.interp1d(x, sig_nn_x_vs_t[-1], kind='linear')
f_sig_yy = interpolate.interp1d(x, sig_yy_x_vs_t[-1], kind='linear')
f_f_y = interpolate.interp1d(x, f_y_x_vs_t[step], kind='linear')  # UNUSED below

# Evaluate σ_nn and σ_yy at x = π/q and x = 2π/q (one minimum, one maximum for cosine)
sig_nn_at_min = f_sig_nn(np.pi / params['q'])
sig_nn_at_max = f_sig_nn(2 * np.pi / params['q'])

sig_yy_at_min = f_sig_yy(np.pi / params['q'])
sig_yy_at_max = f_sig_yy(2 * np.pi / params['q'])

# Store these values back into the parameter dataframe and overwrite CSV
df['sig_nn_at_min'] = [sig_nn_at_min]
df['sig_nn_at_max'] = [sig_nn_at_max]

df['sig_yy_at_min'] = [sig_yy_at_min]
df['sig_yy_at_max'] = [sig_yy_at_max]

df.to_csv('./inputParams_add.csv', index=False)
