import pandas as pd
import pygmmis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#import matplotlib.lines as lines
#import matplotlib.cm
#import datetime
#import logging

df = pd.read_csv("/Users/pfeuffer/Downloads/debugPLFQ/BSAConsensus/consensus_eval/pandas.csv")
df.dropna(inplace=True)
coords = df.as_matrix(['comet_xcorr','msgf_raw'])
gmm = pygmmis.GMM(K=2, D=2)      # K components, D dimensions
logL, U = pygmmis.fit(gmm, coords)  # logL = log-likelihood, U = association of data to components


def plotResults(orig, data, gmm, patch=None, description=None, disp=None):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # plot inner and outer points
    #ax.plot(orig[:,0], orig[:,1], 'o', mfc='None', mec='r', mew=1)
    missing = np.isnan(data)
    if missing.any():
        data_ = data.copy()
        data_[missing] = -25 # put at limits of plotting range
    else:
        data_ = data
    #ax.plot(data_[:,0], data_[:,1], 's', mfc='b', mec='None')#, mew=1)

    # prediction
    B = 100
    x,y = np.meshgrid(np.linspace(0,2.5,B), np.linspace(-100,300,B))
    coords = np.dstack((x.flatten(), y.flatten()))[0]

    # compute sum_k(p_k(x)) for all x
    p = gmm(coords).reshape((B,B))
    # for better visibility use arcshinh stretch
    p = np.arcsinh(p/1e-4)
    cs = ax.contourf(p, 10, extent=(0,2.5,-100,300), cmap=plt.cm.Greys)
    for c in cs.collections:
        c.set_edgecolor(c.get_facecolor())

    # plot boundary
    if patch is not None:
        import copy
        if hasattr(patch, '__iter__'):
            for p in patch:
                ax.add_artist(copy.copy(p))
        else:
            ax.add_artist(copy.copy(patch))

    if description is not None:
        # add description and complete data logL to plot
        logL = gmm(orig, as_log=True).mean()
        ax.text(0.05, 0.95, r'%s' % description, ha='left', va='top', transform=ax.transAxes, fontsize=20)
        ax.text(0.05, 0.89, '$\log{\mathcal{L}} = %.3f$' % logL, ha='left', va='top', transform=ax.transAxes, fontsize=20)

    # show size of error dispersion as Circle
    if disp is not None:

        circ1 = patches.Circle((12.5, -2.5), radius=disp, fc='b', ec='None', alpha=0.5)
        circ2 = patches.Circle((12.5, -2.5), radius=2*disp, fc='b', ec='None', alpha=0.3)
        circ3 = patches.Circle((12.5, -2.5), radius=3*disp, fc='b', ec='None', alpha=0.1)
        ax.add_artist(circ1)
        ax.add_artist(circ2)
        ax.add_artist(circ3)
        ax.text(12.5, -2.5, r'$\sigma$', color='w', fontsize=20, ha='center', va='center')

    ax.set_xlim(0, 2.5)
    ax.set_ylim(-100, 300)
    #ax.set_xticks([])
    #ax.set_yticks([])
    #fig.subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99)
    fig.show()


plotResults(None, coords, gmm)

#test = gmm(coords)
#print(test)