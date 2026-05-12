import pandas as pd
import numpy as np
import csv
from scipy import interpolate
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import vtuIO
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.colors as colors
import proplot as pplt
import os
from copy import deepcopy
from scipy import interpolate
from scipy import integrate
import glob
## SET SOME GENERAL PARAMS AND LOAD COLORS YOU WANT TO USE ###
cc = ['#F3DECA','#FA9483','#2D4057','#4097aa','#453f3d']
cc_special = ['#5eac82']
params = {'subplots.share': False, 'subplots.span': False,'metawidth':1.5,'gridwidth':0,'legend.fontsize':8,'font.sans-serif':'Source Sans Pro'}#'mathtext.fontset':'sourcesanspro'}
params_dark = {'subplots.share': False, 'subplots.span': False,'metawidth':1.5,'gridwidth':0,'legend.fontsize':8,'font.sans-serif':'Source Sans Pro',
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
               "axes.facecolor":"black",
               "figure.facecolor": "black",
               "figure.edgecolor": "black",
               "savefig.facecolor": "black",
               "savefig.edgecolor": "black",
               "boxplot.boxprops.color": "white",
               "boxplot.capprops.color": "white",
               "boxplot.flierprops.color": "white",
               "boxplot.flierprops.markeredgecolor": 'k'"white",
               "boxplot.whiskerprops.color": "white"}


vik = pplt.Colormap('vik')
bilbao = pplt.Colormap('bilbao')


__DATAPATH__ = 'Data/'
__SAVEPATH__ = './'

# read csv i.e. the input parameters for this run
df = pd.read_csv('inputParams_add.csv')
params = df.iloc[0].to_dict()

# result = glob.glob('./Data/stresses*.{}'.format('csv'))
# print(result)
no_csv_files = len(glob.glob('./Data/surface_stresses*.{}'.format('csv')))


iteration_steps = np.arange(0,no_csv_files,1)

sig_yy_x_vs_t = []
sig_nn_x_vs_t = []
f_y_x_vs_t = []
x_vs_t = []
shape_deriv_x_vs_t = []
for step in iteration_steps:
    data = pd.read_csv(__DATAPATH__ +"surface_stresses%04d.csv"%(step))
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

### STANDARD LINE PLOT 
fig = pplt.figure(share=False,aspect=1.61)
# Panel 1
ax1 = fig.subplot(111)
# additional lines
# ax1.axhline(0, linestyle='--', color='black', lw=1)# horizontal line
#legend standard

cycle_reds = pplt.Cycle('reds', no_csv_files,left=0.1)
cycle_blues = pplt.Cycle('blues', no_csv_files,left=0.1)
cycle_greens = pplt.Cycle('greens', no_csv_files,left=0.1)
cc_custom_reds= cycle_reds.by_key()['color']
cc_custom_blues= cycle_blues.by_key()['color']
cc_custom_greens= cycle_greens.by_key()['color']


xnew = np.arange(x_vs_t[0][0], x_vs_t[0][-1]+0.01, 0.01)
for step in iteration_steps:
    x = x_vs_t[step]

    
    f_sig_nn = interpolate.interp1d(x, sig_nn_x_vs_t[step],kind='linear')
    ax1.plot(xnew,f_sig_nn(xnew),color=cc_custom_blues[step],label=r'$\sigma_{nn}$')

f_sig_yy = interpolate.interp1d(x, sig_yy_x_vs_t[-1],kind='linear')
f_f_y = interpolate.interp1d(x, f_y_x_vs_t[-1],kind='linear')
ax1.plot(xnew,f_f_y(xnew),color=cc_custom_greens[-1],label=r'$T_{y}$')
ax1.plot(xnew,f_sig_yy(xnew),color=cc_custom_reds[-1],label=r'$\sigma_{yy}$')


legend_elements = [Line2D([0], [0], color=cc_custom_blues[-1], lw=2, label=r'$\sigma_{nn}$'),
                   Line2D([0], [0], color=cc_custom_greens[-1], lw=2, label=r'$T_{y}$'),
                   Line2D([0], [0], color=cc_custom_reds[-1], lw=2, label=r'$\sigma_{yy}$') ]
ax1.legend(handles=legend_elements,loc='t', ncols=3, frame=False)
# ax1.legend(loc='t', ncols=3, frame=False)
# ax1.spines[['right', 'top']].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.format(xlabel=r'x', ylabel=r'stress (force)',xlim=(x_vs_t[0][0],x_vs_t[0][-1]))#,ylim=(0,1))#
        
fig.save(__SAVEPATH__+'surface_stresses.pdf')
fig.save(__SAVEPATH__+'surface_stresses.png')





### INTEGRAL OF Y-TRACTION
fig = pplt.figure(share=False,aspect=1.61)
# Panel 1
ax1 = fig.subplot(111)
# additional lines
# ax1.axhline(0, linestyle='--', color='black', lw=1)# horizontal line
#legend standard

cycle_grays = pplt.Cycle('grays', no_csv_files,left=0.1)
cc_custom_grays= cycle_grays.by_key()['color']
xnew = np.arange(x_vs_t[0][0], x_vs_t[0][-1]+0.01, 0.01)
for step in iteration_steps:
    x = x_vs_t[step]
    x = x_vs_t[-1]
    f_f_y = interpolate.interp1d(x, f_y_x_vs_t[step],kind='cubic')

    y_int = integrate.cumtrapz(f_y_x_vs_t[step]*np.sqrt(1+shape_deriv_x_vs_t[step]**2), x,initial=0)
    # y_int = integrate.cumtrapz(f_f_y(xnew), xnew,initial=0)

    # ax1.plot(xnew,y_int,color=cc_custom_grays[step],label=r'$\int\;T_y dx$')
    ax1.plot(x,y_int,color=cc_custom_grays[step],label=r'$\int\;T_y dx$')


ax1.axhline(0, linestyle='--', color='black', lw=1)# horizontal line
ax1.axvline(2, linestyle='--', color='red', lw=1)
legend_elements = [Line2D([0], [0], color=cc_custom_grays[-1], lw=2, label=r'$\int\;T_y dx$') ]
ax1.legend(handles=legend_elements,loc='t', ncols=3, frame=False)
# ax1.legend(loc='t', ncols=3, frame=False)
# ax1.spines[['right', 'top']].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.format(xlabel=r'x', ylabel=r'stress (force)',xlim=(x_vs_t[0][0],x_vs_t[0][-1]))#,ylim=(0,1))#
        
fig.save(__SAVEPATH__+'integral_T_y.pdf')
fig.save(__SAVEPATH__+'integral_T_y.png')



# ### save steady state variables
# f_sig_nn = interpolate.interp1d(x, sig_nn_x_vs_t[-1],kind='linear')
# f_sig_yy = interpolate.interp1d(x, sig_yy_x_vs_t[-1],kind='linear')
# f_f_y = interpolate.interp1d(x, f_y_x_vs_t[step],kind='linear')

# last step (steady state)
x = x_vs_t[-1]
xmin, xmax = float(np.min(x)), float(np.max(x))

n = int(params["n"])
L = xmax - xmin
q = 2*np.pi*n / L

# --- ALL CRESTS (cos = +1) ---
k_min_crest = int(np.ceil(q * xmin / (2*np.pi)))
k_max_crest = int(np.floor(q * xmax / (2*np.pi)))
k_crests = np.arange(k_min_crest, k_max_crest + 1)

x_maxpos_all = (2*np.pi * k_crests) / q

# --- ALL TROUGHS (cos = -1) ---
k_min_trough = int(np.ceil((q*xmin - np.pi) / (2*np.pi)))
k_max_trough = int(np.floor((q*xmax - np.pi) / (2*np.pi)))
k_troughs = np.arange(k_min_trough, k_max_trough + 1)

x_minpos_all = (2*np.pi * k_troughs + np.pi) / q

# --- evaluate stresses ---
sig_nn_at_max_all = f_sig_nn(x_maxpos_all)
sig_nn_at_min_all = f_sig_nn(x_minpos_all)

sig_yy_at_max_all = f_sig_yy(x_maxpos_all)
sig_yy_at_min_all = f_sig_yy(x_minpos_all)
# minima (troughs)
df_min = pd.DataFrame({
    'x_pos': x_minpos_all,
    'type': ['min'] * len(x_minpos_all),
    'sig_nn': sig_nn_at_min_all,
    'sig_yy': sig_yy_at_min_all
})

# maxima (crests)
df_max = pd.DataFrame({
    'x_pos': x_maxpos_all,
    'type': ['max'] * len(x_maxpos_all),
    'sig_nn': sig_nn_at_max_all,
    'sig_yy': sig_yy_at_max_all
})

# combine
df_extrema = pd.concat([df_min, df_max], ignore_index=True)

# optional: sort by position
df_extrema = df_extrema.sort_values(by='x_pos')

# save to new CSV
df_extrema.to_csv('./extrema_stresses.csv', index=False)
# df = df.drop(columns=df.columns[0])


# --- mean values over all extrema ---
sig_nn_min_mean = float(np.mean(sig_nn_at_min_all))
sig_nn_max_mean = float(np.mean(sig_nn_at_max_all))

sig_yy_min_mean = float(np.mean(sig_yy_at_min_all))
sig_yy_max_mean = float(np.mean(sig_yy_at_max_all))
df['sig_nn_at_min'] = [sig_nn_min_mean]
df['sig_nn_at_max'] = [sig_nn_max_mean]

df['sig_yy_at_min'] = [sig_yy_min_mean]
df['sig_yy_at_max'] = [sig_yy_max_mean]
df.to_csv('./inputParams_add.csv',index=False)

