import numpy as np
import matplotlib.pyplot as plt


MPROTON_AU = 1836.15
EV_HA = 27.211
ANG_BOHR = .529177
GPA_EV_ANG3 = 160.2176487 # GPa/(eV/Ang^3)

#====================================================================
# SECT 1: library functions
#====================================================================
def ideal_gas( V, kT, xmol, mmol ):

    # Convert length scale from Angstroms to atomic units (bohr)
    logL = 1.0/3*np.log( V/ANG_BOHR**3 )

    # log-canonical partition function
    logQgasmol = 3./2*(2*logL + np.log((kT/EV_HA)*(mmol*MPROTON_AU)/(2*np.pi)))+1.0

    Fgasmol = -kT*logQgasmol;
    Sgasmol = logQgasmol + 3.0/2;

    Fgas = np.sum(xmol*Fgasmol);
    Sgas = np.sum(xmol*Sgasmol);

    return Fgas, Sgas

#====================================================================
def ideal_mix( kT, xHS ):
    Smix_a = -xHS*np.log( xHS )
    Smix_a[xHS==0] = 0.0
    Smix = np.sum( Smix_a )
    Fmix = -kT*Smix
    return Fmix, Smix

#====================================================================
def hard_sphere_mix( V, xHS, dHS, debug_output=False ):
    # hard_sphere_mix( V, xHS, dHS, Soutput=False, debug_output=False ):

    VHS = np.pi/6*dHS**3
    fpackHS = xHS*VHS/V
    fpack = np.sum( fpackHS )

    fpackHSi = np.expand_dims(fpackHS,axis=0)
    dHSi = np.expand_dims(dHS,axis=0)
    xHSi = np.expand_dims(xHS,axis=0)
    # calculate weighted average quantities for ij HS pairs
    dHSij = np.sqrt( dHSi*dHSi.T )

    fpackHSij = np.sqrt(fpackHSi*fpackHSi.T)/fpack
    deltaHSij = fpackHSij * ((dHSi-dHSi.T)/dHSij)**2 * np.sqrt(xHSi*xHSi.T)

    # Note that y1ij and y2ij are symmetric matrices
    y1ij = deltaHSij*(dHSi+dHSi.T)/dHSij
    dHSwtavg = dHSij*np.sum( (fpackHS/fpack)/dHS )
    y2ij = deltaHSij*dHSwtavg

    # sum for all i and j>i
    y1 = 0.5*(np.sum(y1ij)-np.trace(y1ij))
    y2 = 0.5*(np.sum(y2ij)-np.trace(y2ij))
    y3 = np.sum((xHS*(fpackHS/fpack)**2)**(1.0/3))**3

    # calculate excess helmholtz free energy
    FexHS_kT = - 3.0/2*(1-y1+y2+y3)+(3*y2+2*y3)/(1-fpack) \
        + 3.0/2*(1-y1-y2-y3/3)/(1-fpack)**2 \
        + (y3-1)*np.log(1-fpack)

    output = ()
    output += (FexHS_kT,)

    # if (Soutput==True)|(debug_output==True):
    if (debug_output==True):
        # calculate compressibility factor Z
        Z = (1 + fpack + fpack**2 - 3*fpack*(y1+y2*fpack) - y3*fpack**3) \
            *(1-fpack)**(-3)
        S_k = -FexHS_kT + np.log(Z)

        SexHS = S_k - (3-2*fpack)/(1-fpack)**2 + 3 \
            - np.log( (1+fpack+fpack**2-fpack**3)/(1-fpack)**3 )

        debug_data = {'Z':Z, 'S_k':S_k, 'SexHS':SexHS}

    # if (Soutput==True):
    #     output += (SexHS,)

    if (debug_output==True):
        output += (debug_data,)

    return output

#====================================================================
if __name__ == "__main__":
    main()
