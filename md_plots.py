
import matplotlib.pyplot as plt
import numpy as np
import utils
from scipy.stats import gaussian_kde
from IPython import embed
####################################################
#####             Plotting Trajectories        #####
####################################################
'''
Example usage: 
title = "Simulations of {} with {}".format(receptor, ligands)
fname = "ligand-rmsd-alignedtm1to4_{}_{}".format(receptor, ligands.replace(' ', '-'))

ap8_rmsds = [rmsd_from_initial('noh resname AP8', molid) for molid in range(5) if molecule.numframes(molid)]
ax = setup_time_trace(title, 'AP8 RMSD')
add_time_trace(ax, ap8_rmsds, 'AP8')
plt.savefig(fname+'_AP8_trace.png')
plt.ylim(0, 8)
plt.clf()
'''

def _simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

# Trace
def add_time_trace(ax, A, labels = '', colors = 'k', smoothing = 20,
                   alpha=.2, alpha_smooth=1, lw=1, lw_smooth=2):
    if type(colors) == str:
        colors = [colors] * len(A)
    if type(labels) == str:
        labels = [labels] * len(A)
    for dists, color, label in zip(A, colors, labels):
        ax.plot(dists, c = color, lw = lw, alpha = alpha)
        #embed()
        if smoothing:
            ax.plot(utils._sliding_mean(dists, smoothing), c = color, lw = lw_smooth, alpha = alpha_smooth, label = label)

def setup_time_trace(title, ylabel,
                     figsize = (9, 4), axissize=16, titlesize=20, dpi = 300):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_ylabel("{}".format(ylabel), size=axissize) #($\AA$)
    ax.set_xlabel('Time (ns)', size=axissize)
    ax.set_title(title, size=titlesize)
    _simpleaxis(ax)
    return ax

# Kernel Density Estimate   
def _kerneldensityestimate(dataset, limits, n_points):
    if limits is None:
        limits = (dataset.min(), dataset.max())
    x_vals = np.linspace(limits[0], limits[1], n_points)

    if dataset.min() < x_vals.min():
        print("Warning: x values do no include minimum.")
    if dataset.max() > x_vals.max():
        print("Warning: x values do no include maximum.")

    kde = gaussian_kde(dataset)
    return x_vals, kde.evaluate(x_vals)

def add_kde(ax, A, label, burnin = 0, limits=None, n_points=100,
                  color = 'k', lw = 4, flip = False):
    A_eq = np.hstack(dists[burnin:] for dists in A)
    x, y = _kerneldensityestimate(A_eq, limits, n_points)
    if flip: x, y = y, x
    ax.plot(x, y, color, lw=lw)

def add_kde_reps(ax, A, label, burnin = 0, limits=None, n_points=100,
                       colors = 'k', lw = 2,  alpha = 0.3, flip = False):
    A_eq = [dists[burnin:] for dists in A]
    limits = (np.hstack(A_eq).min(), np.hstack(A_eq).max())
    if type(colors) == str:
        colors = [colors] * len(A)
    for dists, color in zip(A_eq, colors):
        x, y = _kerneldensityestimate(dists, limits, n_points)
        if flip: x, y = y, x
        ax.plot(x, y, color = color, lw=lw, alpha = alpha)

def setup_kde(title, xlabel,
                    figsize = (5, 4), axissize=16, titlesize=16, dpi = 300,
                    flip = False):
    ylabel = "Frequency"
    if flip:
        ylabel, xlabel = '', ylabel
        title = ''
        figsize = (figsize[1], figsize[0])
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_ylabel(ylabel, size=axissize)
    ax.set_xlabel(xlabel, size=axissize)
    ax.set_title(title, size=titlesize)
    _simpleaxis(ax)
    return ax
    
# Combined
def setup_trace_and_kde(title, ylabel,
                         figsize = (12, 4), width_ratio = [4, 1], axissize=16, titlesize=20, dpi = 300):
    fig, (trace, kde) = plt.subplots(1, 2, sharey = True, figsize=figsize,
                           gridspec_kw = {'width_ratios':width_ratio}, dpi=dpi)
    trace.set_ylabel("{}".format(ylabel), size=axissize)
    trace.set_xlabel('Time (ns)', size=axissize)
    kde.set_xlabel('Frequency', size=axissize)
    kde.set_xticks([], [])
    trace.set_title(title, size=titlesize)
    _simpleaxis(trace)
    _simpleaxis(kde)
    fig.subplots_adjust(wspace=0)
    return trace, kde