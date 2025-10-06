import scipy.io
import scipy.ndimage
import numpy as np
import matplotlib.pyplot as plt
import sys
import os 
import pandas as pd
import time
import mpl_toolkits.axes_grid1.anchored_artists as mpl
import json


# from sklearn.manifold import Isomap
# from sklearn.decomposition import PCA


os.chdir(r"mpfc_spike_tracking_data") 

files = os.listdir()[:9]


def getnameofactivecells(filename,threshold):
    obj = pd.read_pickle(filename)

    binsize = 0.100   # seconds
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

    return list(np.array(obj['cell_names'])[goodcells]) # returns a list of cells with # of spikes exceeding the threshold

def getsharedcells(threshold): # NOTE:Remember to adjust file numbers for different rats!
    
    sharedcells = set(getnameofactivecells(files[0],threshold))
    for file in files[1:9]:
        sharedcells = sharedcells & set(getnameofactivecells(file,threshold))
    return sharedcells



##### Time-varying % of coactive pairs (Diana chasing 1) ######
obj = pd.read_pickle(files[1])

# Create a matrix of relevant cells and their binned spikes
    
nbins = len(obj['time_bins'])-1
binsize = np.mean(obj['time_bins'][1:]-obj['time_bins'][:(-1)])

goodcellnames = set(np.load('../variables/goodcellnames.npy'))
cellnames = list(getsharedcells(115) & goodcellnames) # Intersection of the cells that are active at least "threshold" number of times in each session and only contain the activity of a single unit. Threshold of 115 results in 313 cells for Diana.

cellindex = [obj['cell_names'].index(f'{cell}') for cell in cellnames]
cellindex.sort()
cellnames = [obj['cell_names'][i] for i in cellindex]
ncells = len(cellnames)

binnedcellspikes = np.zeros((ncells, nbins))
for i in range(ncells):
    cell = obj['cell_activities'][cellindex[i]]
    active_bins = [int(spiketime // binsize) for spiketime in cell]
    for t in active_bins:
        binnedcellspikes[i][t]+=1

######################
# Smoothing
######################
   
sigma = 5 # window for smoothing
for i in range(ncells):
    binnedcellspikes[i,:] = scipy.ndimage.gaussian_filter1d(binnedcellspikes[i,:], sigma)
    binnedcellspikes[i,:] = (binnedcellspikes[i,:]-np.mean(binnedcellspikes[i,:]))/np.std(binnedcellspikes[i,:])

##########################
#Binarization
##########################

# Create a binary coactivity vector for each cell; 1 in bin i if at least 5 bins 
# in the 21 window {i-10,...,i,...,i+10} is greater or equal to 1    ---quantiled value of over 0.90--- 


biactivity = np.zeros(np.shape(binnedcellspikes))
for i in range(ncells):
   for j in range(nbins):
      N=0
      if j>10:
        for k in range(21):
          try:
            if binnedcellspikes[i][j-10+k]>0.90:
              N +=1
            else:
              continue
          except IndexError:
            pass
          continue
        if N>=5:
          biactivity[i][j]=1
        else:
          continue 
      else:
        for k in range(0,j+11):
          if binnedcellspikes[i][k]>0.90:
            N +=1
          else:
            continue
        if N>=5:
          biactivity[i][j]=1
        else:
          continue
np.save('../variables/Diana chasing2 biactivity.npy', biactivity)

biactivity = np.load('../variables/Diana chasing2 biactivity.npy')

####################
#Compute coactive pairs
####################


###### Computes which pairs of cells are significantly coactive ###### Took about 5.5 hrs to run locally. Result saved on local pc.

spikeindices = [np.where(cell==1)[0] for cell in biactivity]
coactivity_score = np.zeros((ncells,ncells))
for i in range(ncells):
   for j in range(ncells):
      rootofproduct = np.sqrt(len(spikeindices[i])*len(spikeindices[j]))
      coactivity_score[i][j] = len(np.intersect1d(spikeindices[i],spikeindices[j]))/rootofproduct - rootofproduct/nbins

num_shuffles = 1000
binary_significance_matrix = np.zeros((ncells,ncells))
 
def get_binary_significance_matrix():
    start=time.time()
    for i in range(ncells):
        for j in range(i,ncells):
            shuffled_coactivity_scores=np.zeros(num_shuffles)
            rootofproduct = np.sqrt(np.sum(biactivity[i])*np.sum(biactivity[j]))
            scaledrootofproduct = rootofproduct/nbins
            for k in range(num_shuffles):
                shift_amount=np.random.randint(1,nbins)
                shuffled_cell=np.roll(biactivity[j],shift_amount)
                shuffled_coactivity_scores[k] = len(np.intersect1d(spikeindices[i],np.where(shuffled_cell>0)[0]))/rootofproduct - scaledrootofproduct
            if np.sum((shuffled_coactivity_scores < coactivity_score[i][j]  )*1) /1000 >= 0.95:
                binary_significance_matrix[i][j] = 1
            else:
                continue
    end=time.time()
    print(end-start)
    np.save('../variables/Diana chasing2 coactive matrix.npy',binary_significance_matrix)
    return

##################


binary_significance_matrix = np.load('../variables/Diana chasing2 coactive matrix.npy')

######## Computes coactivity over time ######## 
start = time.time()
sequence_of_coactivity_matrices=[np.zeros((ncells,ncells),dtype='int8')  for _ in range(nbins)] # takes almost 6 minutes
for i in range(ncells):
  for j in range(i,ncells):
    if binary_significance_matrix[i][j]==1:
       for t in range(nbins):
          if biactivity[i][t]==1 and biactivity[j][t]==1:
             sequence_of_coactivity_matrices[t][i][j]=1
          else:
            continue
    else:
      continue
end = time.time()
duration = end-start

# def get_coactivity_over_time():
#   total_coactive_pairs = np.sum(binary_significance_matrix)  # of unique pairs
#   proportion_of_coactivity_over_time = np.zeros(nbins)
#   for t in range(nbins):
#     proportion_of_coactivity_over_time[t] = np.sum(sequence_of_coactivity_matrices[t])/total_coactive_pairs
#   return proportion_of_coactivity_over_time
  
# coactivity_over_time = get_coactivity_over_time()
# np.save('../variables/Diana chasing2 coactivity over time.npy', coactivity_over_time)
##############

coactivity_over_time = np.load('../variables/Diana chasing2 coactivity over time.npy')



#Plot with x-axis given by bins

# fig, ax = plt.subplots()

# ax.plot(np.arange(nbins), coactivity_over_time, 'k',linewidth=0.3)
# ax.set_xticks([])
# ax.set_xlabel('')
# ax.set_yticks([0,round(max(coactivity_over_time),2)], ['0',f'{max(coactivity_over_time)//0.01}'],fontsize=16)
# ax.set_ylabel('% Co-active neurons',fontsize=16)

# scalebar = mpl.AnchoredSizeBar(ax.transData,nbins/10, '1 min', 'lower left',pad=0.01,color='black',bbox_to_anchor=[250,40.0],frameon=False,size_vertical=0.0001,borderpad=0)
# ax.add_artist(scalebar);
# ax.margins(x=0)
# plt.show()

#Plot with x-axis given by seconds
fig, ax = plt.subplots()

sessionId = ' '.join(obj['output_file_prefix'].split('_')[:2]).capitalize()
sessionId = f"{sessionId[:-1]} {sessionId[-1]}"

ax.plot(np.linspace(0,obj['time_bins'][-1],nbins), coactivity_over_time, 'k',linewidth=0.7)
#ax.set_xticks([''])
ax.set_title(sessionId)
ax.set_xlabel('')
ax.set_yticks([0,round(max(coactivity_over_time),2)], ['0',f'{max(coactivity_over_time)//0.01}'],fontsize=16)
ax.set_ylabel('% Co-active neurons',fontsize=16)

scalebar = mpl.AnchoredSizeBar(ax.transData,60, '1 min', 'lower left',pad=0.01,color='black',bbox_to_anchor=[250,40.0],frameon=False,size_vertical=0.0001,borderpad=0)
ax.add_artist(scalebar)
ax.margins(x=0)
plt.show()