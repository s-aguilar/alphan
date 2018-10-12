import json
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Analyzor.EnergyAndRangeConvertor import PowerTable


def dist(x1,y1,x2,y2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2)


# put it into a pandas data frame
with open('dataV2.dat') as f:
    l = json.load(f)
    info = []
    for d in l:
        x = d['data']
        info.append({'runID': d['runID'], 'fileID': d['fileID'], 'eventID': d['eventID'],
               'x1':x[0][0],'y1':x[0][1],
               'x2':x[1][0],'y2':x[1][1],
               'x3':x[2][0],'y3':x[2][1],
               'xc':x[3][0],'yc':x[3][1]})

    df = pd.DataFrame(info)


# extract the values into list
x =[]
y = []
for index, row in df.iterrows():
    x.append([row["x1"],row["x2"],row["x3"],row["xc"]])
    y.append([row["y1"],row["y2"],row["y3"],row["yc"]])



# define vertex as center and correct the other positions and also calculate ranges (from vertex)
xcentered = []
ycentered = []
ranges = []
posAndRange = []

for index in range(len(x)):
    # center them
    xcentered.append(x[index]-x[index][3])
    ycentered.append(y[index]-y[index][3])
    lengths = dist(xcentered[index].tolist(),ycentered[index].tolist(),xcentered[index][3],ycentered[index][3])
    ranges.append(lengths)

    # note for each event the stuff should be grouped in groups of 3 belonging to each point of triangle
    for step in range(3):
        posAndRange.append([xcentered[index][step],ycentered[index][step],ranges[index][step]])

    # plt.scatter(xcentered[index],ycentered[index],marker ='.',color ='r')
    # plt.scatter(xcentered[index][3],ycentered[index][3],marker='*',color='b')
    # plt.grid(b='True',color='b',alpha=0.8)
    # plt.show()



# lets sort them by lengths of the tracks. Assume the longest track is the initial 10C beam.
#  The second longest track is scattered 10C beam, the shortest track is the scattered
#  4He target. Also convert range to an energy
# tuples are of form:
#   [ [(x1,x2,x3),(y1,y2,y3),(r1,r2,r3),(e1,e2,e3)], [...] ]
sortRange = []
################################# Fuse converting range to energy in this loop############################
for index in range(0,len(posAndRange),3):
    collector = []
    for step in range(3):
        collector.append((posAndRange[index][step],posAndRange[index+1][step],posAndRange[index+2][step]))
    #print(collector)
    #print(sorted(collector, key=lambda rain: rain[2]))
    sortRange.append(sorted(collector, key=lambda rain: rain[2]))   # sort each event by increasing length of tracking.
    print(sortRange)
    break


# Assume the longest track is the initial 10C beam. The second longest track is scattered 10C beam, the shortest track
# is scattered 4He target convert range to an energy
pt10C = PowerTable("Tables/EnergyTable/10C_HeCO2_400_Torr.txt")
pt4He = PowerTable("Tables/EnergyTable/He_HeCO2_400_Torr.txt")

# convert range to energy
posAndEnergy = []
print(sortRange[0])
print(sortRange[0][2])
print(sortRange[0][2][2])

for index in range(len(sortRange)):


    break

#energy = pt10C.R2E()
