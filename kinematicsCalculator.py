import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Analyzor.KINEMATICS.AtomicMassTable import GetElement


# Constants
_c = 3.0e8
_U = 931.4940954 #MeV/c^2
_m4He = GetElement(2,4)[3]*_U
_m10B = GetElement(5,10)[3]*_U
_m11B = GetElement(5,11)[3]*_U
_m10C = GetElement(6,10)[3]*_U


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

def quadraticSolver(a,b,c):
    x1 = (-b+np.sqrt(b*b-4*a*c))/(2*a)
    x2 = (-b-np.sqrt(b*b-4*a*c))/(2*a)
    return x1,x2


def energyToAngleAndEnergy(e_vertex,e1,mass1,ex_energy1=0,ex_energy2=0):
    """ Apply conservation of energy and momentum to find energy and angle """

    # Conserve energy
    deltaE = ex_energy1 + ex_energy2
    e2 = e_vertex - e1 - deltaE

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

    # Apply conservation of momentum
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)

    theta1 = np.arccos((mom0*mom0+mom1*mom1-mom2*mom2)/(2*mom0*mom1)) * 180/np.pi
    theta2 = np.arccos((mom0*mom0+mom2*mom2-mom1*mom1)/(2*mom0*mom2)) * 180/np.pi

    # Check for physical values of angle
    check1 = np.isnan(theta2)

    # Mask values that are not physical
    e1 = e1[check1==0]
    e2 = e2[check1==0]
    theta1 = theta1[check1==0]
    theta2 = theta2[check1==0]

    return e1,theta1,e2,theta2

def angleToAngleEnergy1(e_vertex,theta1,ex_energy1=0,ex_energy2=0):
    """ Given the scattering angle of the alpha particle, apply conservation of
    energy and momentum to find energy and angle
    """

    mass1 = _m4He
    mass2 = _m10C

    # Convert input angle from degrees to radians
    theta1 = theta1 * np.pi/180

    # Find Energy of alpha
    deltaE = ex_energy1 + ex_energy2
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)
    # mom1 = 2*mom0*np.cos(theta1)/(1+mass2/mass1)
    mom1_1,mom1_2 = quadraticSolver(mass2/mass1,-2*mom0*np.cos(theta1),2*mass2*deltaE)
    e1_1 = mom1_1*mom1_1/2/mass1
    e1_2 = mom1_2*mom1_2/2/mass1

    # Find energy of 10C
    e2_1 = e_vertex - e1_1 - deltaE # ex_energy term added
    e2_2 = e_vertex - e1_2 - deltaE # ex_energy term added
    mom2_1 = classicalEnergyToMomentum(e2_1,mass2)
    mom2_2 = classicalEnergyToMomentum(e2_2,mass2)

    # Find other angle
    theta2_1 = np.arcsin(mom1_1*np.sin(theta1)/mom2_1) * 180/np.pi
    theta2_2 = np.arcsin(mom1_2*np.sin(theta1)/mom2_2) * 180/np.pi

    theta1 = np.concatenate((theta1*180/np.pi,theta1*180/np.pi))
    theta2 = np.concatenate((theta2_1,theta2_2))
    e1 = np.concatenate((e1_1,e1_2))
    e2 = np.concatenate((e2_1,e2_2))

    return e1,theta1,e2,theta2

def angleToAngleEnergy2(e_vertex,theta1,ex_energy1=0,ex_energy2=0): ## need to account for inelastic excitation
    """ Given the scattering angle of the alpha particle, apply conservation of
    energy and momentum to find energy and angle
    """

    mass1 = _m4He
    mass2 = _m10C

    theta1 = theta1 * np.pi/180

    # Apply conservation of momentum
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)

    mom1_1,mom1_2 = quadraticSolver(1+mass2/mass1,-2*mom0*np.cos(theta1),ex_energy1+ex_energy2)
    mom2_1 = np.sqrt(2*mass2*(mom0*mom0/2/mass2-mom1_1*mom1_1/2/mass1-ex_energy1-ex_energy2))
    mom2_2 = np.sqrt(2*mass2*(mom0*mom0/2/mass2-mom1_2*mom1_2/2/mass1-ex_energy1-ex_energy2))


    # Assign energies
    e1_1 = mom1_1*mom1_1/2/mass1
    e1_2 = mom1_2*mom1_2/2/mass1
    e2_1 = mom2_1*mom2_1/2/mass1
    e2_2 = mom2_2*mom2_2/2/mass1

    # Find other angle
    theta2_1 = np.arcsin(mom1_1*np.sin(theta1)/mom2_1) * 180/np.pi
    theta2_2 = np.arcsin(mom1_2*np.sin(theta1)/mom2_2) * 180/np.pi

    e1 = np.concatenate((e1_1,e1_2))
    e2 = np.concatenate((e2_1,e2_2))
    theta1 = np.concatenate((theta1,theta1))
    theta2 = np.concatenate((theta2_1,theta2_2))

    return e1,theta1,e2,theta2

def generalAngleToAngleEnergy(e_vertex,mass_vertex,theta_1,mass_1,mass_2):
    """ Given the scattering angle of any particle and its mass and the
    incoming beam energy and species, apply conservation of energy and momentum
    to find energy and angle of outgoing
    """

    theta1 = theta_1 * np.pi/180

    # Apply conservation of momentum
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)

    # 2 solutions from quadratic equation
    mom1_1,mom1_2 = quadraticSolver(1+mass_2/mass_1,-2*mom0*np.cos(theta1),-mom0*mom0*(mass_2/mass_vertex-1))
    mom2_1 = np.sqrt(2*mass_2*(mom0*mom0/2/mass_vertex-mom1_1*mom1_1/2/mass_1))
    mom2_2 = np.sqrt(2*mass_2*(mom0*mom0/2/mass_vertex-mom1_2*mom1_2/2/mass_1))

    # Assign energies
    e1_1 = mom1_1*mom1_1/2/mass_1
    e1_2 = mom1_2*mom1_2/2/mass_1
    e2_1 = mom2_1*mom2_1/2/mass_2
    e2_2 = mom2_2*mom2_2/2/mass_2

    # Find other angle
    theta2_1 = np.arcsin(mom1_1*np.sin(theta1)/mom2_1) * 180/np.pi
    theta2_2 = np.arcsin(mom1_2*np.sin(theta1)/mom2_2) * 180/np.pi

    # Concatenate both solutions from quadratic equation
    e1 = np.concatenate((e1_1,e1_2))
    e2 = np.concatenate((e2_1,e2_2))
    theta1 = np.concatenate((theta1,theta1))
    theta2 = np.concatenate((theta2_1,theta2_2))
    return e1,theta1,e2,theta2

"""
##########
## MAIN ##
##########
"""

# Initialize arrays for 10C,11B,10B
e1_plot,ang1_plot,e2_plot,ang2_plot = np.array([]),np.array([]),np.array([]),np.array([])
_e1_plot,_ang1_plot,_e2_plot,_ang2_plot = np.array([]),np.array([]),np.array([]),np.array([])
__e1_plot,__ang1_plot,__e2_plot,__ang2_plot = np.array([]),np.array([]),np.array([]),np.array([])

erange = np.linspace(34,35,3)

for iii in erange:
    # reset input angle array after each iteration
    ang1, _ang1, __ang1 = np.linspace(0,90,91),np.linspace(0,90,91),np.linspace(0,90,91)

    ev = iii*np.ones(91)

    # 10C
    e1,ang1,e2,ang2 = generalAngleToAngleEnergy(ev,_m10C,ang1,_m4He,_m10C)
    e1_plot = np.append(e1_plot,e1)
    ang1_plot = np.append(ang1_plot,ang1)
    e2_plot = np.append(e2_plot,e2)
    ang2_plot = np.append(ang2_plot,ang2)

    # 11B
    _e1,_ang1,_e2,_ang2 = generalAngleToAngleEnergy(ev,_m11B,_ang1,_m4He,_m11B)
    _e1_plot = np.append(_e1_plot,_e1)
    _ang1_plot = np.append(_ang1_plot,_ang1)
    _e2_plot = np.append(_e2_plot,_e2)
    _ang2_plot = np.append(_ang2_plot,_ang2)

    # 10B
    __e1,__ang1,__e2,__ang2 = generalAngleToAngleEnergy(ev,_m10B,__ang1,_m4He,_m10B)
    __e1_plot = np.append(__e1_plot,__e1)
    __ang1_plot = np.append(__ang1_plot,__ang1)
    __e2_plot = np.append(__e2_plot,__e2)
    __ang2_plot = np.append(__ang2_plot,__ang2)


# 3D surface plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(e1_plot,e2_plot,ang2_plot,color='b', alpha=0.8)#,label=r'$^{10}C(\alpha,\alpha)^{10}C$')
# ax.plot(_e1_plot,_e2_plot,_ang2_plot,color='r', alpha=0.8,label=r'$^{11}B(\alpha,\alpha)^{11}B$')
# ax.plot(__e1_plot,__e2_plot,__ang2_plot,color='y', alpha=0.8,label=r'$^{10}B(\alpha,\alpha)^{10}B$')

ax.set_xlabel("e1")
ax.set_ylabel("e2")
ax.set_zlabel(r"$\theta_{2}$")
ax.legend()
plt.show()


# plt.scatter(_ang1,_ang2)
# plt.title(r"$\theta_{2}$ vs $\theta_{1}$")
# plt.xlabel(r"$\theta_{1} [deg]$")
# plt.ylabel(r"$\theta_{2} [deg]$")
# plt.show()
#
# plt.scatter(_ang1,_e1)
# plt.title(r"$E_{1}$ vs $\theta_{1}$")
# plt.xlabel(r"$\theta_{1} [deg]$")
# plt.ylabel("E$_{1}$ [MeV]")
# plt.show()
#
# plt.scatter(_ang2,_e2)
# plt.title(r"$E_{2}$ vs $\theta_{2}$")
# plt.xlabel(r"$\theta_{2} [deg]$")
# plt.ylabel("E$_{2}$ [MeV]")
# plt.show()
#
# plt.scatter(e1,e2)
# # plt.ylim(0,15)
# # plt.xlim(0,15)
# plt.show()
