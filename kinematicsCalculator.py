import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from Analyzor.KINEMATICS.tbjcconstants import *


# Constants
_c = 3.0e8
_U = 931.4940954 #MeV/c^2
_m4He = 4.00260305404663*_U
_m10C = 10.0168533325195*_U

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
    e2 = e_vertex - e1 - deltaE # ex_energy term added

    # Sort
    if mass1 == _m4He:
        mass2 = _m10C
    elif mass1 == _m10C:
        mass2 = _m4He
    else:
        raise Exception("It is neither a helium nor a carbon")

    # Assign momentums
    mom1 = classicalEnergyToMomentum(e1,mass1) # ex_energy term added
    mom2 = classicalEnergyToMomentum(e2,mass2) # ex_energy term added

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

    massb = mass_vertex # this one changes

    mass1 = mass_1
    mass2 = mass_2

    # if mass_1 == _m4He:
    #     mass1 = mass_1
    #     mass2 = _m10C
    # elif mass_1 == _m10C:
    #     mass1 = mass_1
    #     mass2 = _m4He

    theta1 = theta_1 * np.pi/180

    # Apply conservation of momentum
    mom0 = classicalEnergyToMomentum(e_vertex,_m10C)

    # 2 solutions from quadratic equation
    mom1_1,mom1_2 = quadraticSolver(1+mass2/mass1,-2*mom0*np.cos(theta1),-mom0*mom0*(mass2/massb-1))
    mom2_1 = np.sqrt(2*mass2*(mom0*mom0/2/massb-mom1_1*mom1_1/2/mass1))
    mom2_2 = np.sqrt(2*mass2*(mom0*mom0/2/massb-mom1_2*mom1_2/2/mass1))

    # Assign energies
    e1_1 = mom1_1*mom1_1/2/mass1
    e1_2 = mom1_2*mom1_2/2/mass1
    e2_1 = mom2_1*mom2_1/2/mass2
    e2_2 = mom2_2*mom2_2/2/mass2

    # Find other angle
    theta2_1 = np.arcsin(mom1_1*np.sin(theta1)/mom2_1) * 180/np.pi
    theta2_2 = np.arcsin(mom1_2*np.sin(theta1)/mom2_2) * 180/np.pi

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
#ev = 15*np.ones(91)
#e1 = np.linspace(0,15,91)
ang1 = np.linspace(0,90,91)

# e1,ang1,e2,ang2 = energyToAngleAndEnergy(ev,e1,_m4He,ex_energy1=0,ex_energy2=0)
# e1,ang1,e2,ang2 = angleToAngleEnergy1(ev,ang1,0,1)
#e1,ang1,e2,ang2 = angleToAngleEnergy2(ev,ang1,_m10C)

e1_plot,ang1_plot,e2_plot,ang2_plot = [],[],[],[]
_e1_plot,_ang1_plot,_e2_plot,_ang2_plot = [],[],[],[]

erange = np.linspace(34,35,100)

for iii in erange:
    # rest input angle array after each iteration
    ang1 = np.linspace(0,90,91)
    _ang1 = np.linspace(0,90,91)
    ev = (1+iii)*np.ones(91)
    e1,ang1,e2,ang2 = generalAngleToAngleEnergy(ev,_m10C,ang1,_m4He,_m10C)
    # this second plot will be used to see if curves are separate
    _e1,_ang1,_e2,_ang2 = generalAngleToAngleEnergy(ev,_m10C,_ang1,_m4He,_m10C)
    e1_plot.append(e1)
    ang1_plot.append(ang1)
    e2_plot.append(e2)
    ang2_plot.append(ang2)
    _e1_plot.append(_e1)
    _ang1_plot.append(_ang1)
    _e2_plot.append(_e2)
    _ang2_plot.append(_ang2)

# 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(e1_plot,e2_plot,ang2_plot)
ax.scatter(_e1_plot,_e2_plot,_ang2_plot)

ax.set_xlabel('e1')
ax.set_ylabel('e2')
ax.set_zlabel(r'$\ttheta_{2}')
plt.show()


# plt.scatter(ang1,ang2)
# plt.title(r"$\theta_{2}$ vs $\theta_{1}$")
# plt.xlabel(r"$\theta_{1} [deg]$")
# plt.ylabel(r"$\theta_{2} [deg]$")
# plt.show()
#
# plt.scatter(ang1,e1)
# plt.title(r"$E_{1}$ vs $\theta_{1}$")
# plt.xlabel(r"$\theta_{1} [deg]$")
# plt.ylabel("E$_{1}$ [MeV]")
# plt.show()
#
# plt.scatter(ang2,e2)
# plt.title(r"$E_{2}$ vs $\theta_{2}$")
# plt.xlabel(r"$\theta_{2} [deg]$")
# plt.ylabel("E$_{2}$ [MeV]")
# plt.show()
#
# plt.scatter(e1,e2)
# # plt.ylim(0,15)
# # plt.xlim(0,15)
# plt.show()
