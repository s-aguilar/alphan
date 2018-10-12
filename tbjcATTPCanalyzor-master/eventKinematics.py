import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from Analyzor.EnergyAndRangeConvertor import PowerTable
from Analyzor.KINEMATICS.KINEMATICS import KINEMATICS
from Analyzor.KINEMATICS.tbjcconstants import *


# Constants
_c = 3.0e8
_U = 931.4940954 #MeV/c^2
_m4He = 4.00260305404663*_U
_m10C = 10.0168533325195*_U


def dist(x1,y1,x2,y2):
    """ Returns the distance between two points (in 2D) """
    return np.sqrt((x1-x2)**2+(y1-y2)**2)
def kineticEnergyAfterCollision(E10):
    """ The projectile is the 10C, target is 4He at rest. Using only
    conservation of kinetic energy, return the kinetic energy of both assuming
    head-on collisions
    """
    ke10C = E10 * ((_m10C-_m4He)/(_m10C+_m4He))**2
    ke4He = E10 * (1 - ((_m10C-_m4He)/(_m10C+_m4He))**2 )
    return [ke10C,ke4He]
def angleBetween(Ax,Ay,Bx,By,Cx,Cy):
    """ Law of Cosines Cos A = (-a^2 + b^2 + c^2)/(2bc). The angle A belongs to
    the vertex, B belongs to beam, C belongs to scattered particle. Return the
    supplementary angle to A
    """
    a = dist(Bx,By,Cx,Cy)
    b = dist(Ax,Ay,Cx,Cy)
    c = dist(Ax,Ay,Bx,By)
    angleInRad = np.arccos( (-a*a+ b*b + c*c)/(2*b*c+0.0001) )
    angle = angleInRad * 180/np.pi
    return 180. - angle
def angleEnergyFinder(e_vertex,e1,mass1):
    """ Apply conservation of energy and momentum to find energy and angle """

    # Conserve energy
    e2 = e_vertex - e1

    # Sort
    if mass1 == _m4He:
        mass2 = _m10C
    elif mass1 == _m10C:
        mass2 = _m4He
    else:
        raise Exception("It is neither a helium nor a carbon")

    # Assign momentums
    mom1 = classicalEnergyToMomentum(e1,mass1)
    mom2 = classicalEnergyToMomentum(e2,mass2)

    #theta1 = theta1*np.pi/180. # convert to radians

    # Apply conservation of momentum
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)

    theta1 = np.arccos((mom0*mom0+mom1*mom1-mom2*mom2)/(2*mom0*mom1)) * 180/np.pi
    theta2 = np.arccos((mom0*mom0+mom2*mom2-mom1*mom1)/(2*mom0*mom2)) * 180/np.pi

    return theta1,e2,theta2
def classicalEnergyToMomentum(energy,mass):
    """ unit: energy(MeV), mass(MeV/c^2)"""
    p2 = 2*mass*energy
    return np.sqrt(p2)
def relativisticEnergyToMomentum(energy,mass):
    """ unit: energy(MeV), mass(MeV/c^2)"""
    p2 = energy*energy - mass*mass
    return np.sqrt(p2)
def kineticEnergyToRelativisticEnergy(KE,mass):
    return KE + mass


###############################################################################
###############################################################################


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
    # WRONG
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
pt10C = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/10C_HeCO2_400_Torr.txt")
pt4He = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/He_HeCO2_400_Torr.txt")

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

##Maybe issue with these guys since some are larger than e10
#print(e1Loss_val)
#print(e2Loss_val)

# e10 is the energy available at the vertex
e_beam = 34. # MeV
ebeam = []
e10 = []
fwhm = 1 # MeV
sigma = fwhm/2.355
for x in range(len(evLoss_val)):
    energy = sigma*np.random.randn()+e_beam
    ebeam.append(energy)
    e10.append(energy-ebeamLoss_val[index])



e10C,e4He = [kineticEnergyAfterCollision(x)[0] for x in e10], \
            [kineticEnergyAfterCollision(x)[1] for x in e10]


# Clean up the data
data = []
for index in range(len(angle1)):
    event = dict(Range=rbeam_val[index],L1=r1_val[index],L2=r2_val[index], \
                 phi1=angle1[index],phi2=angle2[index],e_10=e10[index], \
                 e1_loss=e1Loss_val[index],e2_loss=e2Loss_val[index])
    data.append(event)

df1 = pd.DataFrame(data)

mask = ( (df1['Range']>2) & df1['L1'].between(4,10000) & df1['L2'].between(4,10000)
& df1['phi1'].between(5,160) & df1['phi2'].between(5,160) )
#& (df1['e_10']>0.0) & (df1['e1_loss']+df1['e2_loss']<df1['e_10']) & (df1['phi1']+df1['phi2'] < 95) )

df1_clean = df1[mask]

ran = df1_clean['Range']
l1 = df1_clean['L1']
l2 = df1_clean['L2']
xx = df1_clean['phi1']
yy = df1_clean['phi2']





# plt.hist2d(xx,yy,
# #xx,yy = np.sort(df1_clean[['phi1','phi2']],axis=1).T
# #plt.hist2d(np.concatenate([xx,yy]),np.concatenate([yy,xx]),
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
#
#
#
# # 0 = C10
# # 1 = He4
# # 2 = He4
# # 3 = C10
# KIN = KINEMATICS((10,4,4,10),40,0,0)
# angles2,angles3 = [],[]
#
# for i in range(181):
#     KIN.calculate(i*constant_degree,0)
#     angles2.append(KIN.thetalab2)
#     angles3.append(KIN.thetalab3)
#
# plt.scatter(angles2,angles3, color='b',s=5)
# plt.scatter(angles3,angles2, color='b',s=5)
# plt.show()
