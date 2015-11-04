import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
from scipy import interpolate

fact = scipy.misc.factorial

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
def hard_sphere_PDF( V, xHS, dHS, rmax=5.0, N=101 ):
    """
    Calculate contact value of pair distribution function
    """
    dHSmax = np.max( dHS )
    r = np.linspace( 0, rmax*dHSmax, N )
    dHSij = 0.5*(dHS+np.expand_dims(dHS,1))
    gij_contact = hard_sphere_contact_PDF( V, xHS, dHS )

    lapfun0 = lambda s, i, j, V=V, xHS=xHS, dHS=dHS:\
        hard_sphere_LT_PDF( s, V, xHS, dHS )[i,j]

    pdf = np.zeros((len(xHS),len(xHS),len(r)))
    for i in np.arange(len(dHS)):
        for j in np.arange(len(dHS)):
            lapfun0ij = lambda s, i=i, j=j: lapfun0(s,i,j)
            rfiltij = r[r>dHSij[i,j]]
            ynuminvij_filt = inv_laplace_euler( lapfun0ij,
                                               rfiltij-dHSij[i,j],
                                               tshft=dHSij[i,j] )
            gij_filt = ynuminvij_filt/rfiltij
            ind_real = np.where( ~np.isnan(gij_filt) )[0]
            gij_filt, rfiltij = ( np.hstack((gij_contact[i,j],gij_filt[ind_real])),
                                np.hstack((dHSij[i,j],rfiltij[ind_real])) )

            gfunc = interpolate.interp1d( rfiltij, gij_filt, kind='cubic' )
            gij = np.zeros( r.shape )
            gij[r>=dHSij[i,j]] = gfunc(r[r>=dHSij[i,j]])
            pdf[i,j,:] = gij

    return pdf, r

#====================================================================
def hard_sphere_contact_PDF( V, xHS, dHS ):
    rho = 1.0/V
    dHSij = 0.5*(dHS+np.expand_dims(dHS,1))
    rhoi = xHS*rho
    fpack = np.pi/6*np.sum(rhoi*dHS**3)

    zeta = fpack*xHS*dHS**3
    gii_contact = ((1-zeta) + 1.5*fpack*np.sum(xHS*dHS**2)*dHS)/(1-zeta)**2

    dgii = np.expand_dims(dHS,1)*gii_contact
    gij_contact = 0.5*(dgii+dgii.T)/dHSij
    return gij_contact

#====================================================================
def hard_sphere_LT_PDF( s, V, xHS, dHS ):
    """
    Calculate Laplace Transform of Pair Distribution Function from Percus-Yevick

     - s must be a scalar currently
     - calculate for i,j HS pair
     - Percus-Yevick approx retrieved by setting alpha=L2=0 from Yuste(1998)
    """
    nHS = len(xHS)

    rho = 1.0/V
    dHSij = 0.5*(dHS+np.expand_dims(dHS,1))
    rhoi = xHS*rho
    fpack = np.pi/6*np.sum(rhoi*dHS**3)

    lam0 = 2*np.pi/(1-fpack)
    lam1 = np.pi**2 * np.sum(rhoi*dHS**2)/(1-fpack)**2

    # Calc intermediate matrices
    # These are much simpler in std. PY vs full Yuste formalism
    L0 = np.tile(lam0 + lam1*dHS, (nHS,1))
    L1 = lam0*dHSij + 0.5*lam1*dHS*np.expand_dims(dHS,1)

    L = L0 + L1*s


    def phi1(x):
        return ((1.0-x)-np.exp(-x))/x**2
    def phi2(x):
        return ((1.0-x+0.5*x**2)-np.exp(-x))/x**3


    A = np.tile(rhoi,(nHS,1)).T \
        *( L0 * np.tile(phi2(dHS*s)*dHS**3,(nHS,1)).T
          + L1 * np.tile(phi1(dHS*s)*dHS**2,(nHS,1)).T )

    Aadj = np.eye(A.shape[0])-A
    XSOLVE = np.linalg.solve(Aadj.T, L.T).T
    # NOTE L2ij=0 yields Percus-Yevick approx
    G = np.exp(-s*dHSij)/(2*np.pi*s**2)*XSOLVE
    return G

#====================================================================
def inv_laplace_euler( lapfun0, t, tshft=1, Nterms=60, meuler=30):
    discreteErrOrd = 10
    A = discreteErrOrd*np.log(10)

    lapfun = lambda s: np.exp(s*tshft)*lapfun0(s)

    # tminArr = tshft*np.logspace(-4,-1,30)
    # funCheck = lapfun(A./2./tminArr + 1j*np.pi)
    # tmin = tminArr[~np.isnan(funCheck)][0]
    # t = linspace(tmin,tmax-tshft,Nsamp);

    tMat = np.tile(t,(Nterms,1))
    kMat = np.tile(np.arange(0,Nterms),(len(t),1)).T

    scomplex = (A+2*kMat*np.pi*1j)/(2*tMat)
    shp = scomplex.shape
    lapMat = np.zeros(scomplex.shape)
    for ind,scomplexi in enumerate(scomplex.ravel()):
        ij = np.unravel_index( ind, shp)
        lapMat[ij[0],ij[1]] = np.squeeze(np.real(lapfun(scomplexi)))

    snMat = (-1)**kMat/tMat*lapMat

    snMat[0,:] = 0.5*snMat[0,:]
    snCum = np.exp(A/2)*np.cumsum(snMat,axis=0)

    ind=np.arange(0,meuler+1)
    bincoeff=fact(meuler)/fact(ind)/fact(meuler-ind)*2**(-meuler)
    wtMat = np.tile(bincoeff,(len(t),1)).T

    EsumAvg = np.zeros(t.shape)
    Nmax = snCum.shape[0]
    numterm = np.arange(meuler,Nmax-meuler+1)
    for i in np.arange(0,len(numterm)):
        inumterm=numterm[i]
        Esum = np.sum(wtMat*snCum[inumterm-1:inumterm+meuler+1,:],axis=0)
        EsumAvg = i/(i+1)*EsumAvg + 1./(i+1)*Esum;

    tVal = t+tshft;
    invLapVal = EsumAvg;
    return invLapVal

#====================================================================
def csteh(n, i):
    acc = 0.0
    for k in xrange(int(np.floor((i+1)/2.0)), int(min(i, n/2.0))+1):
        num = k**(n/2.0) * fact(2 * k)
        den = fact(i - k) * fact(k -1) * fact(k) * fact(2*k - i) * fact(n/2.0 - k)
        acc += (num /den)
    expo = i+n/2.0
    term = np.power(-1+0.0j,expo)
    res = term * acc
    return res.real

#====================================================================
def nlinvsteh(F, t, n = 6):
    acc = 0.0
    lton2 = np.log(2) / t
    for i in xrange(1, n+1):
        a = csteh(n, i)
        b = F(i * lton2)
        acc += (a * b)
    return lton2 * acc
#====================================================================
if __name__ == "__main__":
    main()
