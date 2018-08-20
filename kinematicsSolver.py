import numpy as np

# Constants
_c = 3.0e8
_U = 931.4940954 #MeV/c^2
_m4He = 4.00260305404663*_U
_m10C = 10.0168533325195*_U

def dist(x1,y1,x2,y2):
    """ Returns the distance between two points (in 2D) """
    return np.sqrt((x1-x2)**2+(y1-y2)**2)

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

def kineticEnergyAfterCollision(E10):
    """ The projectile is the 10C, target is 4He at rest. Using only
    conservation of kinetic energy, return the kinetic energy of both assuming
    head-on collisions
    """
    ke10C = E10 * ((_m10C-_m4He)/(_m10C+_m4He))**2
    ke4He = E10 * (1 - ((_m10C-_m4He)/(_m10C+_m4He))**2 )
    return [ke10C,ke4He]

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

    # Concatenate both solutions from quadratic equation but suppress the solutions that are 0
    if mom1_1.all() == 0.:
        #print("suppress zeroes1")
        return e1_2,theta1,e2_2,theta2_2
    elif mom1_2.all() ==0.:
        #print("supress zeroes2")
        return e1_1,theta1,e2_1,theta2_1
    else:
        #print("no zeroes")
        e1 = np.concatenate((e1_1,e1_2))
        e2 = np.concatenate((e2_1,e2_2))
        theta1 = np.concatenate((theta1,theta1))
        theta2 = np.concatenate((theta2_1,theta2_2))
        return e1,theta1,e2,theta2
