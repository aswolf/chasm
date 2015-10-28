import numpy as np
import matplotlib.pyplot as plt


MPROTON_AU = 1836.15
EV_HA = 27.211
ANG_BOHR = .529177
GPA_EV_ANG3 = 160.2176487 # GPa/(eV/Ang^3)


#====================================================================
# SECT 1: library functions
#====================================================================
def calc_closest_packing_dist( V, xHS, NoxyHS ):
    Noxy = xHS*NoxyHS
    Voxy = V/Noxy
    fpack_hcp = np.pi/np.sqrt(18)
    distoxy = (6/np.pi*fpack_hcp*Voxy)**(1.0/3)
    return distoxy

def calc_ideal_packing( V, xHS, NoxyHS ):
    Noxy = xHS*NoxyHS
    Voxy = V/Noxy
    fpack_hcp = np.pi/np.sqrt(18)
    distoxy = (6/np.pi*fpack_hcp*Voxy)**(1.0/3)
    return distoxy

