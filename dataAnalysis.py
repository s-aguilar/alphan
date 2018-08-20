import json
import numpy as np
import pandas as pd
from kinematicsSolver import *
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from Analyzor.KINEMATICS.tbjcconstants import *
from Analyzor.KINEMATICS.KINEMATICS import KINEMATICS
from Analyzor.EnergyAndRangeConvertor import PowerTable
from Analyzor.KINEMATICS.AtomicMassTable import GetElement

# Constants
_c = 3.0e8
_U = 931.4940954 #MeV/c^2
_m4He = GetElement(2,4)[3]*_U
_m10B = GetElement(5,10)[3]*_U
_m11B = GetElement(5,11)[3]*_U
_m10C = GetElement(6,10)[3]*_U


"""
##########
## MAIN ##
##########
"""
###############
# THEORETICAL #
###############

# Initialize arrays for 10C,11B,10B
e1_plot,ang1_plot,e2_plot,ang2_plot,r1_plot,r2_plot = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
_e1_plot,_ang1_plot,_e2_plot,_ang2_plot,_r1_plot,_r2_plot = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
__e1_plot,__ang1_plot,__e2_plot,__ang2_plot,__r1_plot,__r2_plot = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

erange = np.linspace(8,15,8)

pt10C = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/10C_HeCO2_400_Torr.txt")
pt4He = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/He_HeCO2_400_Torr.txt")

for iii in erange:
    # reset input angle array after each iteration
    ang1, _ang1, __ang1 = np.linspace(0,90,91),np.linspace(0,90,91),np.linspace(0,90,91)

    ev = iii*np.ones(91)

    # 10C
    e1,ang1,e2,ang2 = generalAngleToAngleEnergy(ev,_m10C,ang1,_m4He,_m10C)
    e1_plot = np.append(e1_plot,e1)
    r1_plot = np.append(r1_plot,pt4He.E2R(e1))
    ang1_plot = np.append(ang1_plot,ang1)
    e2_plot = np.append(e2_plot,e2)
    r2_plot = np.append(r2_plot,pt10C.E2R(e2))
    ang2_plot = np.append(ang2_plot,ang2)

    # 11B
    _e1,_ang1,_e2,_ang2 = generalAngleToAngleEnergy(ev,_m11B,_ang1,_m4He,_m11B)
    _e1_plot = np.append(_e1_plot,_e1)
    _r1_plot = np.append(_r1_plot,pt4He.E2R(_e1))
    _ang1_plot = np.append(_ang1_plot,_ang1)
    _e2_plot = np.append(_e2_plot,_e2)
    _r2_plot = np.append(_r2_plot,pt10C.E2R(_e2))
    _ang2_plot = np.append(_ang2_plot,_ang2)

    # 10B
    __e1,__ang1,__e2,__ang2 = generalAngleToAngleEnergy(ev,_m10B,__ang1,_m4He,_m10B)
    __e1_plot = np.append(__e1_plot,__e1)
    __r1_plot = np.append(__r1_plot,pt4He.E2R(__e1))
    __ang1_plot = np.append(__ang1_plot,__ang1)
    __e2_plot = np.append(__e2_plot,__e2)
    __r2_plot = np.append(__r2_plot,pt10C.E2R(__e2))
    __ang2_plot = np.append(__ang2_plot,__ang2)


################
# EXPERIMENTAL #
################

# Load data into a pandas data frame
with open('/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/data/full_data.dat') as f:
    l = json.load(f)
    info = []
    for d in l:
        x = d['data']
        info.append({'runID':d['runID'],'fileID':d['fileID'],'eventID':d['eventID'],
               'x1':x[0][0],'y1':x[0][1],
               'x2':x[1][0],'y2':x[1][1],
               'x3':x[2][0],'y3':x[2][1],
               'xc':x[3][0],'yc':x[3][1]})

    df = pd.DataFrame(info)



# extract the values into list, calculate range and append
x =[]
ranges = []
y = []
for index, row in df.iterrows():
    r1 = dist(row["x1"],row["y1"],row["xc"],row["yc"])
    r2 = dist(row["x2"],row["y2"],row["xc"],row["yc"])
    r3 = dist(row["x3"],row["y3"],row["xc"],row["yc"])
    rc = np.float64(0.)
    x.append([row["x1"],row["x2"],row["x3"],row["xc"]])
    y.append([row["y1"],row["y2"],row["y3"],row["yc"]])
    ranges.append([r1,r2,r3,rc])
    #WRONG
    # r1 = dist(row["x1"],row["y1"],row["xc"],row["yc"])
    # r2 = dist(row["x2"],row["y3"],row["xc"],row["yc"])
    # r3 = dist(row["x2"],row["y3"],row["xc"],row["yc"])
    # rc = np.float64(0.)
    # x.append([row["x1"],row["x2"],row["x3"],row["xc"]])
    # y.append([row["y1"],row["y2"],row["y3"],row["yc"]])
    # ranges.append([r1,r2,r3,rc])


# Assume the right most point belongs to the incoming beam 10C. Sorting makes
# the right most point become the last tuple in the list for each event.
#
# Tuples are of form:
#   [ [(x1,y1,r1),(x2,y2,r2),(x3,y3,r3),(x4,y4,r4),], [...] ]
# pt10C = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/10C_HeCO2_400_Torr.txt")
# pt4He = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/He_HeCO2_400_Torr.txt")

kinematicInfo = []

for index in range(len(x)):
    posAndRange = []

    # Note for each event the stuff should be grouped in groups of 4 belonging
    # to each point of triangle and vertex
    for step in range(4):
        posAndRange.append((x[index][step],y[index][step],ranges[index][step]))

    # First sort by x value in increasing order (last element is beam point),
    # Next sort the 3 smaller x values by range so vertex is first
    sortedRangeX = sorted(posAndRange,key=lambda z: z[0])      # sort by x
    sortedRangeR = sorted(sortedRangeX[:-1],key=lambda z: z[2]) # sort by range
    sortedRangeR.append(sortedRangeX[3])
    sortedFinal = sortedRangeR

    # Convert range to energy and add to tuples
    eVert = 0.
    eloss10C= float(pt10C.R2E(sortedFinal[1][2]))
    eloss4He = float(pt4He.R2E(sortedFinal[2][2]))
    elossBeam = float(pt10C.R2E(sortedFinal[3][2]))

    sortedFinal[0] = sortedFinal[0] + (eVert,)
    sortedFinal[1] = sortedFinal[1] + (eloss10C,)
    sortedFinal[2] = sortedFinal[2] + (eloss4He,)
    sortedFinal[3] = sortedFinal[3] + (elossBeam,)

    # Calculate angle and append
    theta1 = angleBetween(sortedFinal[0][0],sortedFinal[0][1], \
                          sortedFinal[3][0],sortedFinal[3][1], \
                          sortedFinal[1][0],sortedFinal[1][1])
    theta2 = angleBetween(sortedFinal[0][0],sortedFinal[0][1], \
                          sortedFinal[3][0],sortedFinal[3][1], \
                          sortedFinal[2][0],sortedFinal[2][1])

    sortedFinal[1] = sortedFinal[1] + (theta1,)
    sortedFinal[2] = sortedFinal[2] + (theta2,)

    kinematicInfo.append(sortedFinal)



# Tuples are of form:
#   [ [(x0,y0,r0,e0),(x1,y1,r1,e1,t1),(x2,y2,r2,e2,t2),(x3,y3,r3,e3)], [...] ]
#
# Everything is sorted now unpack them into lists and plot them
xv_val = [x[0][0] for x in kinematicInfo]
x1_val = [x[1][0] for x in kinematicInfo]
x2_val = [x[2][0] for x in kinematicInfo]
xbeam_val = [x[3][0] for x in kinematicInfo]

yv_val = [x[0][1] for x in kinematicInfo]
y1_val = [x[1][1] for x in kinematicInfo]
y2_val = [x[2][1] for x in kinematicInfo]
ybeam_val = [x[3][1] for x in kinematicInfo]

rv_val = [x[0][2] for x in kinematicInfo]
r1_val = [x[1][2] for x in kinematicInfo]
r2_val = [x[2][2] for x in kinematicInfo]
rbeam_val = [x[3][2] for x in kinematicInfo]

evLoss_val = [x[0][3] for x in kinematicInfo]
e1Loss_val = [x[1][3] for x in kinematicInfo]
e2Loss_val = [x[2][3] for x in kinematicInfo]
ebeamLoss_val = [x[3][3] for x in kinematicInfo]

angle1 = [x[1][4] for x in kinematicInfo]
angle2 = [x[2][4] for x in kinematicInfo]

# e10 is the energy available at the vertex
e_beam = 34.5 # MeV
ebeam = []
e10 = []
fwhm = 1 # MeV
sigma = fwhm/2.355
for x in range(len(evLoss_val)):
    energy = sigma*np.random.randn()+e_beam
    ebeam.append(energy)
    e10.append(energy-ebeamLoss_val[x])


e10C,e4He = [kineticEnergyAfterCollision(x)[0] for x in e10], \
            [kineticEnergyAfterCollision(x)[1] for x in e10]
####################################################################
# plt.scatter(pt10C.E2R(e10C),pt4He.E2R(e4He))
plt.scatter(e10C,e4He)
plt.show()
#####################################################################
# Clean up the data
data = []
for index in range(len(angle1)):
    event = dict(Range=rbeam_val[index],L1=r1_val[index],L2=r2_val[index], \
                 phi1=angle1[index],phi2=angle2[index],e_10=e10[index], \
                 e1_loss=e1Loss_val[index],e2_loss=e2Loss_val[index])
    data.append(event)

df1 = pd.DataFrame(data)

mask = ( (df1['Range']>2) & df1['L1'].between(4,10000) & df1['L2'].between(4,10000)
      & df1['phi1'].between(5,160) & df1['phi2'].between(5,160) & (df1['e_10']>0) )

df1_clean = df1[mask]

ran = df1_clean['Range']
l1 = df1_clean['L1']
e1_loss = df1_clean['e1_loss']
e2_loss = df1_clean['e2_loss']
l2 = df1_clean['L2']
phi1 = df1_clean['phi1']
phi2 = df1_clean['phi2']
e_10 = df1_clean['e_10']


# plt.hist2d(phi1,phi2,
# # phi1,phi2 = np.sort(df1_clean[['phi1','phi2']],axis=1).T
# # plt.hist2d(np.concatenate([phi1,phi2]),np.concatenate([phi2,phi1]),
#            bins=[np.arange(0,100,1),np.arange(0,100,1)],
#            norm=LogNorm());
# poly = [[60,35],[35,60],[30,50],[50,30],[60,35]]
# path2 = mplPath.Path(poly)
# patch2 = patches.PathPatch(path2,facecolor=(0,0,0,0),EdgeColor='r',lw=3,\
#                            linestyle='--')
# plt.gca().add_patch(patch2)
#
#
# poly = [[75,5],[5,75],[5,5],[75,5]]
# path1 = mplPath.Path(poly)
# patch1 = patches.PathPatch(path1,facecolor=(0,0,0,0),EdgeColor='b',lw=3,\
#                            linestyle='--')
# plt.gca().add_patch(patch1)
# plt.show()




# 3D surface plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(r1_plot,r2_plot,ang2_plot,color='b', alpha=0.8)#,label=r'$^{10}C(\alpha,\alpha)^{10}C$')
ax.scatter(l2,l1,phi1,color='r', s=0.5, alpha = 0.5)
ax.set_xlim(0,50)
ax.set_ylim(0,50)
ax.set_zlim(0,90)
# ax.plot(_e1_plot,_e2_plot,_ang2_plot,color='r', alpha=0.8,label=r'$^{11}B(\alpha,\alpha)^{11}B$')
# ax.plot(__e1_plot,__e2_plot,__ang2_plot,color='y', alpha=0.8,label=r'$^{10}B(\alpha,\alpha)^{10}B$')

ax.set_xlabel("range 1")
ax.set_ylabel("range 2")
ax.set_zlabel(r"$\theta_{2}$")
ax.legend()
plt.show()
