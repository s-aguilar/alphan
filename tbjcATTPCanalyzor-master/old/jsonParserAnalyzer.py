import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import seaborn as sns
import matplotlib.patches as patches
from Analyzor import Func
from Analyzor.EnergyAndRangeConvertor import PowerTable
from Analyzor.KINEMATICS.KINEMATICS import KINEMATICS
from Analyzor.KINEMATICS.tbjcconstants import *


plt.rcParams['figure.figsize'] = (20, 20)


with open('data.dat') as f:
    l = json.load(f)
    data = []
    for one in l:
        print(one)
        # print(one['data'][0][2])
        # break
        xc,yc = one['data'][0][2]
        points = one['data'][0][:-1]
        d = [Func.GetLineInfo((x_,y_),(xc,yc)) for x_,y_ in points]
        data.append({'runID': one['runID'], 'fileID': one['fileID'], 'eventID': one['eventID'],
               'theta1':d[0][0],'L1':d[0][1],
               'theta2':d[1][0],'L2':d[1][1],'Range':d[2][1]})

df = pd.DataFrame(data)

mask = ((df['Range']>2) & df['L1'].between(4,10000) & df['L2'].between(4,10000)
& df['theta1'].between(5,160) & df['theta2'].between(5,160))

df_tmp = df[mask]

x,y = np.sort(df_tmp[['theta1','theta2']],axis=1).T
#print np.concatenate([x,y])
plt.hist2d(np.concatenate([x,y]),np.concatenate([y,x]),
           bins=[np.arange(0,100,1),np.arange(0,100,1)],
           norm=LogNorm());


poly = [[60,35],[35,60],[30,50],[50,30],[60,35]]
path2 = mplPath.Path(poly)
patch2 = patches.PathPatch(path2, facecolor=(0,0,0,0),EdgeColor='r', lw=5,linestyle='--')
plt.gca().add_patch(patch2)


poly = [[75,5],[5,75],[5,5],[75,5]]
path1 = mplPath.Path(poly)
patch1 = patches.PathPatch(path1, facecolor=(0,0,0,0),EdgeColor='b', lw=5,linestyle='--')
plt.gca().add_patch(patch1)







'''

# read data
with open('data.dat') as f:
    data = json.load(f)

# lists (theta,range)
t1 = []
r1 = []
t2 = []
r2 = []
t3 = []
r3 = []
trash = []

# loop through all data
for d in data[:]:
    # get data list
    x = d['data']

    if(x[0][0][0]<180.0 and x[0][1][0]<180.0 and x[0][2][0]<180.0):
        t1.append(x[0][0][0])
        r1.append(x[0][1][0])
        t2.append(x[0][1][0])
        r3.append(x[0][1][1])
        t3.append(x[0][2][0])
        r3.append(x[0][2][1])

pt = PowerTable("Tables/EnergyTable/He_HeCO2_400_Torr.txt")
blah = pt.E2R
blah = []
blah2 = []
print (r1)
blah = PowerTable.GetE0(r1)
print (blah)

'''
'''
    # sort and append to new lists
    if(x[0][0][0]<90.0 and x[0][1][0]<90.0 and x[0][2][0]<90.0):
        continue
    if(x[0][1][0]<90.0 and x[0][0][0]<90.0):
        t3.append(x[0][2][0])
        r3.append(x[0][2][1])
        t1.append(x[0][0][0])
        r1.append(x[0][0][1])
        t2.append(x[0][1][0])
        r2.append(x[0][1][1])
    elif(x[0][2][0]<90.0 and x[0][1][0]<90.0):
        t3.append(x[0][0][0])
        r3.append(x[0][0][1])
        t1.append(x[0][1][0])
        r1.append(x[0][1][1])
        t2.append(x[0][2][0])
        r2.append(x[0][2][1])
    elif(x[0][2][0]<90.0 and x[0][0][0]<90.0):
        t3.append(x[0][1][0])
        r3.append(x[0][1][1])
        t1.append(x[0][2][0])
        r1.append(x[0][2][1])
        t2.append(x[0][0][0])
        r2.append(x[0][0][1])
'''
# cast lists into arrays
'''
t1 = np.array(t1)
r1 = np.array(r1)
t2 = np.array(t2)
r2 = np.array(r2)
t3 = np.array(t3)
r3 = np.array(r3)
'''
# sns.jointplot(t1,t2,kind="kde")
# plt.xlabel('Angle1')
# plt.ylabel('Angle2')
# plt.grid(True)
# plt.show()

'''
plt.hist2d(t1,t3,bins=50,norm=LogNorm())
plt.colorbar()
plt.xlabel('Angle1')
plt.ylabel('Angle2')
plt.grid(True)
plt.show()
'''
'''
plt.hist2d(r1,r2,bins=50,norm=LogNorm())
plt.colorbar()
# plt.plot(r1,r2,linestyle='',marker='o',markersize=0.7)
plt.xlabel('Range1')
plt.ylabel('Range2')
plt.grid(True)
plt.show()

plt.hist2d(t1,r1,bins=50,norm=LogNorm())
plt.colorbar()
plt.xlabel('Theta1')
plt.ylabel('Range1')
plt.grid(True)
plt.show()

plt.hist2d(t2,r2,bins=50,norm=LogNorm())
plt.colorbar()
plt.xlabel('Theta2')
plt.ylabel('Range2')
plt.grid(True)
plt.show()
'''
    # sort list
    # if (x[0][0][0] == 180.0):
    #     x[0].append(x[0][0])
    #     del x[0][0]
    # elif (x[0][1][0] == 180.0):
    #     x[0].append(x[0][1])
    #     del x[0][1]
