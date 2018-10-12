import numpy as np
import matplotlib.pyplot as plt
from Analyzor.KINEMATICS.tbjcconstants import *
from Analyzor.KINEMATICS.KINEMATICS import KINEMATICS
from Analyzor.EnergyAndRangeConvertor import PowerTable



# 0 = 8Li
# 1 = He4
# 2 = n
# 3 = 11B
_m8Li = 8.02248764038086
_m4He = 4.00260305404663
_mn = 1.00866496562958
_m11B = 11.0093050003052

pt11B = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/B_HeCO2_400Torr.txt")

'''
Change the energy inputs in following variables
'''
nrg1 = 1.4
nrg2 = 1.3
nrg3 = 1.5


###############################################################################
###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,0)
r1,e1,ang1 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang1.append(KIN.thetalab3)
    e1.append(KIN.K3)
    r1.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,0)
r11,e11,ang11 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang11.append(KIN.thetalab3)
    e11.append(KIN.K3)
    r11.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,0)
r12,e12,ang12 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang12.append(KIN.thetalab3)
    e12.append(KIN.K3)
    r12.append(pt11B.E2R(KIN.K3))

###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,2.124)
r2,e2,ang2 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang2.append(KIN.thetalab3)
    e2.append(KIN.K3)
    r2.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,2.124)
r21,e21,ang21 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang21.append(KIN.thetalab3)
    e21.append(KIN.K3)
    r21.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,2.124)
r22,e22,ang22 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang22.append(KIN.thetalab3)
    e22.append(KIN.K3)
    r22.append(pt11B.E2R(KIN.K3))

###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,4.445)
r3,e3,ang3 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang3.append(KIN.thetalab3)
    e3.append(KIN.K3)
    r3.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,4.445)
r31,e31,ang31 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang31.append(KIN.thetalab3)
    e31.append(KIN.K3)
    r31.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,4.445)
r32,e32,ang32 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang32.append(KIN.thetalab3)
    e32.append(KIN.K3)
    r32.append(pt11B.E2R(KIN.K3))

###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,5.02)
r4,e4,ang4 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang4.append(KIN.thetalab3)
    e4.append(KIN.K3)
    r4.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,5.02)
r41,e41,ang41 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang41.append(KIN.thetalab3)
    e41.append(KIN.K3)
    r41.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,5.02)
r42,e42,ang42 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang42.append(KIN.thetalab3)
    e42.append(KIN.K3)
    r42.append(pt11B.E2R(KIN.K3))

###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,6.742)
r5,e5,ang5 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang5.append(KIN.thetalab3)
    e5.append(KIN.K3)
    r5.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,6.742)
r51,e51,ang51 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang51.append(KIN.thetalab3)
    e51.append(KIN.K3)
    r51.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,6.742)
r52,e52,ang52 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang52.append(KIN.thetalab3)
    e52.append(KIN.K3)
    r52.append(pt11B.E2R(KIN.K3))

###############################################################################

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg1,0,6.792)
r6,e6,ang6 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang6.append(KIN.thetalab3)
    e6.append(KIN.K3)
    r6.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg3,0,6.792)
r61,e61,ang61 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang61.append(KIN.thetalab3)
    e61.append(KIN.K3)
    r61.append(pt11B.E2R(KIN.K3))

KIN = KINEMATICS((_m8Li,_m4He,_mn,_m11B),nrg2,0,6.792)
r62,e62,ang62 = [],[],[]
for i in range(181):
    KIN.calculate(i*constant_degree,0)
    ang62.append(KIN.thetalab3)
    e62.append(KIN.K3)
    r62.append(pt11B.E2R(KIN.K3))

###############################################################################


# Plot Energy vs Angle
plt.plot(ang1,e1, color='r',label='$E_{x}(^{11}B)=0\ MeV$',linewidth=3)#,linestyle='--')
plt.plot(ang11,e11, color='r',linestyle='--')
plt.plot(ang12,e12, color='r',linestyle='--')

plt.plot(ang2,e2, color='b',label='$E_{x}(^{11}B)=2.124\ MeV$',linewidth=3)
plt.plot(ang21,e21, color='b',linestyle='--')
plt.plot(ang22,e22, color='b',linestyle='--')

plt.plot(ang3,e3, color='m',label='$E_{x}(^{11}B)=4.445\ MeV$',linewidth=3)
plt.plot(ang31,e31, color='m',linestyle='--')
plt.plot(ang32,e32, color='m',linestyle='--')

plt.plot(ang4,e4, color='c',label='$E_{x}(^{11}B)=5.02\ MeV$',linewidth=3)
plt.plot(ang41,e41, color='c',linestyle='--')
plt.plot(ang42,e42, color='c',linestyle='--')

plt.plot(ang5,e5, color='y',label='$E_{x}(^{11}B)=6.742\ MeV$',linewidth=3)
plt.plot(ang51,e51, color='y',linestyle='--')
plt.plot(ang52,e52, color='y',linestyle='--')

plt.plot(ang6,e6, color='k',label='$E_{x}(^{11}B)=6.792\ MeV$',linewidth=3)
plt.plot(ang61,e61, color='k',linestyle='--')
plt.plot(ang62,e62, color='k',linestyle='--')

plt.xlim(0,70)
plt.ylim(0,3.5)
plt.legend()
plt.title('$^{8}Li(\\alpha,n)^{11}B\ \ \ E_{k}(^{8}Li): %6.2f,%6.2f,%6.2f\ MeV$' % (nrg2,nrg1,nrg3))
plt.xlabel('$\Theta_{4}$ [deg]')
plt.ylabel('$E_{4}$ [MeV]')
plt.show()



# Plot Range vs Angle
plt.plot(ang1,r1, color='r',label='$E_{x}(^{11}B)=0\ MeV$',linewidth=3)
plt.plot(ang11,r11, color='r',linestyle='--')
plt.plot(ang12,r12, color='r',linestyle='--')

plt.plot(ang2,r2, color='b',label='$E_{x}(^{11}B)=2.124\ MeV$',linewidth=3)
plt.plot(ang21,r21, color='b',linestyle='--')
plt.plot(ang22,r22, color='b',linestyle='--')

plt.plot(ang3,r3, color='m',label='$E_{x}(^{11}B)=4.445\ MeV$',linewidth=3)
plt.plot(ang31,r31, color='m',linestyle='--')
plt.plot(ang32,r32, color='m',linestyle='--')

plt.plot(ang4,r4, color='c',label='$E_{x}(^{11}B)=5.02\ MeV$',linewidth=3)
plt.plot(ang41,r41, color='c',linestyle='--')
plt.plot(ang42,r42, color='c',linestyle='--')

plt.plot(ang5,r5, color='y',label='$E_{x}(^{11}B)=6.742\ MeV$',linewidth=3)
plt.plot(ang51,r51, color='y',linestyle='--')
plt.plot(ang52,r52, color='y',linestyle='--')

plt.plot(ang6,r6, color='k',label='$E_{x}(^{11}B)=6.792\ MeV$',linewidth=3)
plt.plot(ang61,r61, color='k',linestyle='--')
plt.plot(ang62,r62, color='k',linestyle='--')

plt.xlim(0,80)
plt.ylim(0,5)
plt.legend()
plt.title('$^{8}Li(\\alpha,n)^{11}B\ \ \ E_{k}(^{8}Li): %6.2f,%6.2f,%6.2f\ MeV$' % (nrg2,nrg1,nrg3))
plt.xlabel('$\Theta_{4}$ [deg]')
plt.ylabel('Range [cm]')
plt.show()






# # Test SRIM
# pt11B = PowerTable("/Users/sebastian/Desktop/ActiveTargetGroup/tbjcATTPCanalyzor-master/Tables/EnergyTable/B_HeCO2_400Torr.txt")
# x = pt11B.E2R(3) # input in MeV
# print(x) # output in mm
