import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import seaborn as sns
import matplotlib.path as mplPath
import matplotlib.patches as patches
from Analyzor.EnergyAndRangeConvertor import PowerTable
import math



with open('data.dat') as f:
    l = json.load(f)
    info = []
    for d in l:
        x = d['data']
        info.append({'runID': d['runID'], 'fileID': d['fileID'], 'eventID': d['eventID'],
               'theta1':x[0][0][0],'L1':x[0][0][1],
               'theta2':x[0][1][0],'L2':x[0][1][1],'Range':x[1]})

df = pd.DataFrame(info)

mask = ((df['Range']>2) & df['L1'].between(4,10000) & df['L2'].between(4,10000)
& df['theta1'].between(5,160) & df['theta2'].between(5,160))

df_tmp = df[mask]


print(df_tmp.head())

for index, row in df_tmp.iterrows():
    rad = []
    r = []
    x = 0
    for x in range(4):
        rad.append()


# x,y = np.sort(df_tmp[['theta1','theta2']],axis=1).T
# plt.hist2d(np.concatenate([x,y]),np.concatenate([y,x]),
#            bins=[np.arange(0,100,1),np.arange(0,100,1)],
#            norm=LogNorm());
#
#
# poly = [[60,35],[35,60],[30,50],[50,30],[60,35]]
# path2 = mplPath.Path(poly)
# patch2 = patches.PathPatch(path2, facecolor=(0,0,0,0),EdgeColor='r', lw=5,linestyle='--')
# plt.gca().add_patch(patch2)
#
#
# poly = [[75,5],[5,75],[5,5],[75,5]]
# path1 = mplPath.Path(poly)
# patch1 = patches.PathPatch(path1, facecolor=(0,0,0,0),EdgeColor='b', lw=5,linestyle='--')
# plt.gca().add_patch(patch1)
#
# plt.show()



pt = PowerTable("Tables/EnergyTable/10C_HeCO2_400_Torr.txt")

Energy2 = []
Energy3 = []
# convert range
Energy1 = pt.R2E(df['L1'])
# print(df['L1'])
# print(Energy1)
# print(pt.E2R(Energy1))

#rold = [24.585051358884648,131.37088451531258,41.866833739669175]
#told = [148.07716780862702,17.72697061311414,180.0]






r1 = np.array([67.47496,38.753840000000004,71.60608,69.04872])
t1 = np.array([15.700000000000001,10.4,14.4,15.200000000000001])
r2 = np.array([23.6,5.0,7.9,19.900000000000002])
t2 = np.array([38.16368,47.40952,48.196400000000004,40.52432])
r3 = np.array([15.100000000000001,16.2,9.6,14.9])
t3 = np.array([86.55680000000001,9.04912,38.16368,57.638960000000004])

radians = np.deg2rad(t3)

rad = [radians[0],radians[1],radians[2],radians[3]]
R = [r3[0],r3[1],r3[2],r3[3]]


# ax = plt.subplot(111, projection='polar')
# ax.scatter(rad, R)
# ax.set_rmax(25)
# plt.show()

# x = np.cos(radians)*r3
# y = np.sin(radians)*r3
#
# plt.scatter(x,y)
# plt.show()




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

    # sort and append to new lists
    if(x[1] > 2 and 4 <= x[0][0][1] <= 10000 and 4 <= x[0][1][1] <= 10000 and 5 <= x[0][0][0] <= 160 and 5 <= x[0][1][0] <= 160):
        t1.append(x[0][0][0])
        r1.append(x[0][0][1])
        t2.append(x[0][1][0])
        r2.append(x[0][1][1])
        t3.append(x[0][2][0])
        r3.append(x[0][2][1])
    # if(x[0][0][0]<90.0 and x[0][1][0]<90.0 and x[0][2][0]<90.0):
    #     continue

# cast lists into arrays
t1 = np.array(t1)
r1 = np.array(r1)
t2 = np.array(t2)
r2 = np.array(r2)
t3 = np.array(t3)
r3 = np.array(r3)

# sns.jointplot(t1,t2,kind="kde")
# plt.xlabel('Angle1')
# plt.ylabel('Angle2')
# plt.grid(True)
# plt.show()

x,y = np.sort([t1,t2], axis = 1)

plt.hist2d(np.concatenate([x,y]),np.concatenate([y,x]),
           bins=[np.arange(0,100,5),np.arange(0,100,5)],
           norm=LogNorm());
plt.colorbar()
plt.grid(True)
plt.show()
'''

'''
plt.hist2d(t1,t2,bins=50,norm=LogNorm())
plt.colorbar()
plt.xlabel('Angle1')
plt.ylabel('Angle2')
plt.grid(True)
plt.show()

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
