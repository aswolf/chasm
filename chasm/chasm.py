import numpy as np
import matplotlib.pyplot as plt


MPROTON_AU = 1836.15
EV_HA = 27.211
ANG_BOHR = .529177
GPA_EV_ANG3 = 160.2176487 # GPa/(eV/Ang^3)

#====================================================================
# SECT 1: main script
#====================================================================
def ideal_gas( V, kT, xHS, mHS ):

    # Convert length scale from Angstroms to atomic units (bohr)
    logL = 1.0/3*np.log( V/ANG_BOHR**3 )

    logZgasHS = 3./2*(2*logL + np.log((kT/EV_HA)*(mHS*MPROTON_AU)/(2*np.pi)))+1.0

    FgasHS = -kT*logZgasHS;
    SgasHS = logZgasHS + 3.0/2;

    Fgas = np.sum(xHS*FgasHS);
    Sgas = np.sum(xHS*SgasHS);

    return Fgas, Sgas


#====================================================================
if __name__ == "__main__":
    main()
