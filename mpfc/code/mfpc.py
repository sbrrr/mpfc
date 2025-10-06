import pandas as pd
import sys
import os
import numpy as np

os.chdir(r"mpfc_spike_tracking_data") 
files = os.listdir()


cellcount = []
filename = []

for file in files:
    obj = pd.read_pickle(f"{file}")
    cellcount.append(len(obj['cell_names']))
    filename.append(f"{file}")
dict = {"cellcount": cellcount, "filename": filename}

print(dict)
df = pd.DataFrame(dict)
df

concatednatedcells = [[],[],[]]
for i in range(3):
    for j in range(9):
        obj = pd.read_pickle(files[i*9 + j])

        concatednatedcells[i] += obj['cell_names']

for l in concatednatedcells:
    l = list(dict.fromkeys(l))
    print(len(l))

obj=pd.read_pickle(files[0])
obj.keys()

thresholds = [[],[],[]]
sharedcells = []

# obj = pd.read_pickle(files[18])

for threshold in range(1025,1026):
    obj = pd.read_pickle(files[0])
    target = 313
    sharedcells = []

    for k in range(len(obj['cell_names'])):
        if len(obj['cell_activities'][k]) >= threshold:
            sharedcells.append(obj['cell_names'][k])
    for i in range(1,9):
        obj = pd.read_pickle(files[i])
        goodcells = []
        for j in range(len(obj['cell_names'])):
            if len(obj['cell_activities'][j]) >= threshold:
                goodcells.append(obj['cell_names'][j])
        sharedcells = list(set(sharedcells) & set(goodcells))
        print(len(goodcells))
        print(len(sharedcells), "\n")

    print(len(sharedcells))

    if len(sharedcells) == target:
        thresholds[0].append(threshold)


    # if len(sharedcells) <= target:
    #     break
print(thresholds)

#     if len(obj['cell_activities']) >= 1000:
#         print(obj['cell_names'][i])
# sharedcells = [pd.read_pickle(files[9*i])['cell_names'] for i in range(3)]
# print([len(sharedcells[i]) for i in range(3)])

# for i in range(3):
#     for j in range(1,9):
#         obj = pd.read_pickle(files[i*9 + j])
#         sharedcells[i] = list(set(sharedcells[i]) & set(obj['cell_names']))

# print([len(cells) for cells in sharedcells])

#cheks if sharedcells[i] is contained in the intersection of cells that were active in every session for rat i
for i in range(3):
    for j in range(9):
        set(sharedcells[i]).issubset(set(pd.read_pickle(files[i*9+j])['cell_names']))
    
print(len(sharedcells))

for file in files[:9]:
    set(sharedcells).issubset(set(pd.read_pickle(file)['cell_names']))

    