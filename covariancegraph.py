import pandas as pd
import sys
import os
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.anchored_artists as mpl

os.chdir(r"mpfc_spike_tracking_data")
files = os.listdir()

#Functions for getting the cells 
def getnameofactivecells(filename,threshold):  ### With 100ms bins
    obj = pd.read_pickle(filename)

    binsize = 0.100     # seconds
    nbins = int(obj['time_bins'][-1] // binsize +1) 
    binnedspikes = np.zeros((len(obj['cell_names']), nbins))

    for i in range(len(obj['cell_names'])):
        cell = obj['cell_activities'][i]
        active_bins = [int(time // binsize) for time in cell]
        for j in active_bins:
            binnedspikes[i][j]+=1

    goodcells = np.zeros(len(binnedspikes[:,0]))<1
    for i in range(len(goodcells)):
        if np.sum(binnedspikes[i])<threshold:
            goodcells[i]=False

    return list(np.array(obj['cell_names'])[goodcells]) # returns a list with the name of cells where # of active cells is greater than or equal to threshold

def getsharedcells(threshold): # NOTE:Remember to adjust file numbers for different rats!  [0:9] for Diana, [9:18] for Leia, [18:27] for Rey 
    
    sharedcells = set(getnameofactivecells(files[0],threshold))
    for file in files[1:9]:
        sharedcells = sharedcells & set(getnameofactivecells(file,threshold))
    return sharedcells

obj = pd.read_pickle(files[1])

# Create a matrix of relevant cells and their binned spikes
    
binsize = 0.100
nbins = int(obj['time_bins'][-1] // binsize +1)

goodcellnames = set(np.load('../goodcellnames.npy'))
cellnames = list(getsharedcells(1) & goodcellnames) # Intersection of the cells that are active at least "threshold" number of times in each session and only contain the activity of a single unit. Threshold of 1025 results in 313 cells for Diana.

cellindex = [obj['cell_names'].index(f'{cell}') for cell in cellnames]
cellindex.sort()
cellnames = [obj['cell_names'][i] for i in cellindex]
ncells = len(cellnames)

binnedspikes = np.zeros((ncells, nbins))
for i in range(ncells):
    cell = obj['cell_activities'][cellindex[i]]
    active_bins = [int(spiketime // binsize) for spiketime in cell]
    for t in active_bins:
        binnedspikes[i][t]+=1

# Square rooting spike count
binnedspikes = np.sqrt(binnedspikes)

#Normalization
binnedspikes = np.array([scipy.stats.zscore(binnedspikes[i]) for i in range(ncells)])

#Moving window covariance
windowsize = 10
nwindows = nbins - windowsize + 1
covscore = np.zeros(nwindows)
for t in range(nwindows):
    window = binnedspikes[:,t:t+windowsize]
    #window = window*(1 * window>=0) #Removes negative cov values
    covscore[t] = np.sum(np.cov(np.abs(window)))

#Plotting
sessionId = ' '.join(obj['output_file_prefix'].split('_')[:2]).capitalize()
sessionId = f"{sessionId[:-1]} {sessionId[-1]}"

#plt.plot(np.linspace(0,obj['time_bins'][-1],nwindows),covscore)
fig, ax = plt.subplots()

ax.plot(np.linspace(0,obj['time_bins'][-1],nwindows),covscore, 'k')
#ax.set_xticks([])
ax.set_title(sessionId)
ax.set_xlabel('seconds')
ax.set_yticks([0,round(max(covscore),2)], ['0',f'{round(max(covscore))}'],fontsize=16)
ax.set_ylabel('Covariance',fontsize=16)

scalebar = mpl.AnchoredSizeBar(ax.transData,60, '1 min', 'lower left',pad=0.01,color='black',bbox_to_anchor=[250,40.0],frameon=False,size_vertical=0.0001,borderpad=0)
ax.add_artist(scalebar);
ax.margins(x=0)
plt.show()
