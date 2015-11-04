import numpy as np
import matplotlib.pyplot as plt
import pdb
import hsmix
import scipy as sp

#====================================================================
def test_ideal_gas_press():
    TOL = .03
    xHS = 1.0

    # atmospheric conditions
    mHS = 28.97
    kT = 1.0/40 # 300 K

    V0 = 39270.0
    V_a = V0*np.array([0.99,1.01])

    # 1bar = kT/V0*1.6e6
    # V_a = V0*np.linspace( .5, 5, 1001)

    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
    Fgas_a = np.zeros( V_a.shape )
    Sgas_a = np.zeros( V_a.shape )
    for ind, V in enumerate( V_a ):
        iFgas, iSgas = hsmix.ideal_gas( V, kT, xHS, mHS )
        Fgas_a[ind] = iFgas
        Sgas_a[ind] = iSgas

    P = -np.diff(Fgas_a)/np.diff(V_a)*hsmix.GPA_EV_ANG3

    assert np.abs(np.log(P*1e4/1.013)) < TOL, \
        'Must recover 1 bar atmospheric pressure'

#====================================================================
def test_ideal_gas_entropy():
    TOL = 1e-3
    xHS = 1.0

    # atmospheric conditions
    mHS = 28.97
    kT0 = 1.0/40 # 300 K

    kT_a = kT0*np.array([.99,1.01])
    V = 39270.0

    # 1bar = kT/V0*1.6e6
    # V_a = V0*np.linspace( .5, 5, 1001)

    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
    Fgas_a = np.zeros( kT_a.shape )
    Sgas_a = np.zeros( kT_a.shape )
    for ind, kT in enumerate( kT_a ):
        iFgas, iSgas = hsmix.ideal_gas( V, kT, xHS, mHS )
        Fgas_a[ind] = iFgas
        Sgas_a[ind] = iSgas

    S = -np.diff(Fgas_a)/np.diff(kT_a)

    assert np.abs( np.log( np.mean(Sgas_a)/S ) ) < TOL

    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

#====================================================================
def test_ideal_mix():
    kT = 1.0
    xHS = np.array([.5,.5])
    Fmix, Smix = hsmix.ideal_mix( kT, xHS )

    assert Smix == np.log(2), 'Smix of 50/50 mix should equal log(2)'

    Fmix, Smix = hsmix.ideal_mix( kT, np.array([0.0,1.0]) )
    assert Smix==0, 'Purely 1 component yields Smix=0'

#====================================================================
def test_hard_sphere_mix():
    TOL = 1e-2

    fpackHS_a=np.array([0.2333, 0.2692, 0.3106, 0.3583, 0.3808, 0.4393, 0.5068])
    dHS = np.array([1, 3])
    xHS = np.array([0.5, 0.5])

    V_a = np.sum( xHS*np.pi/6*dHS**3 )/fpackHS_a

    FexHS_kT = np.zeros( V_a.shape )
    debug_output = None

    debug_output = None
    for ind, V in enumerate(V_a):
        iFexHS_kT, idebug_output = hsmix.hard_sphere_mix( V, xHS, dHS, debug_output=True )
        FexHS_kT[ind] = iFexHS_kT
        if debug_output is None:
            debug_output = {}
            for key in idebug_output:
                debug_output[key] = np.array(idebug_output[key])

        else:
            for key in idebug_output:
                debug_output[key] = np.append(debug_output[key],
                                              idebug_output[key])


    Z_a = np.array([2.368,2.772,3.356,4.241,4.764,6.567,9.898])
    Sk_a = -np.array([0.139,.205,.306,.467,.564,.898,1.495])

    assert np.all(np.abs(np.log(debug_output['S_k']/Sk_a)) < TOL),  \
        'S_k  values disagree with Mansoori 1971 Table 2.'
    assert np.all(np.abs(np.log(debug_output['Z']/Z_a)) < TOL), \
        'Z  values disagree with Mansoori 1971 Table 2.'

    assert False,  'excess S values do not match Mansoori 1971 Table 2 values'
    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

#====================================================================
def test_bessel_inv_laplace_euler():
    TOL = 1e-6

    t=np.linspace(1e-3,15,100)

    # Bessel function test (ringing oscillation)
    lapfun0 = lambda s: 1.0/np.sqrt(s**2+1)
    ynuminv = hsmix.inv_laplace_euler( lapfun0, t, tshft=0.0 )
    yexact = sp.special.jv(0,t)

    assert np.all( np.abs(ynuminv-yexact) < TOL ), \
        'numerical inverse not within tolerance'

#====================================================================
def test_invsqrt_cos_inv_laplace_euler():
    TOL = 1e-6

    t=np.linspace(0.1,20,100)

    # Bessel function test (ringing oscillation)
    lapfun0 = lambda s: 1.0/np.sqrt(s)*np.exp(-1.0/s)
    ynuminv = hsmix.inv_laplace_euler( lapfun0, t, tshft=0.0 )
    yexact = 1.0/np.sqrt(np.pi*t)*np.cos(np.sqrt(4*t))

    assert np.all( np.abs(ynuminv-yexact) < TOL ), \
        'numerical inverse not within tolerance'

#====================================================================
def test_exp_cos_inv_laplace_euler():
    TOL = 1e-6

    t=np.linspace(0.1,20,100)
    omega = 2.0
    a = 1.0

    # Bessel function test (ringing oscillation)
    lapfun0 = lambda s, a=a, omega=omega: (s+a)/((s+a)**2+omega**2)

    ynuminv = hsmix.inv_laplace_euler( lapfun0, t, tshft=0.0 )
    yexact = np.exp(-a*t)*np.cos(omega*t)

    assert np.all( np.abs(ynuminv-yexact) < TOL ), \
        'numerical inverse not within tolerance'

#====================================================================
def test_shifted_exp_cos_inv_laplace_euler():
    TOL = 1e-6

    N = 1001
    N = 101

    omega = 4.0
    a = 0.8
    tshft = 1.5

    delt = np.linspace(0.01,6,N)
    t = delt+ tshft
    yexact = np.exp(-a*(t-tshft))*np.cos(omega*(t-tshft))+1.0
    yexact = np.exp(-a*(t-tshft))*np.cos(omega*(t-tshft))+1.0
    yexact[t<tshft] = 0.0

    lapfun0 = lambda s, a=a, omega=omega: \
         np.exp(-s*tshft)*( (s+a)/((s+a)**2+omega**2) + 1.0/s )

    # ynuminv = hsmix.inv_laplace_euler( lapfun0, t, tshft=0.0 )
    ynuminv = hsmix.inv_laplace_euler( lapfun0, delt, tshft=tshft )

    # NOTE nan value at first value
    dely = ynuminv-yexact
    dely = dely[~np.isnan(dely)]

    #plt.clf()
    #plt.plot(t,yexact,'r-',t,ynuminv,'k-')

    assert np.all( np.abs(dely) < TOL ), \
        'numerical inverse not within tolerance'

#====================================================================
def test_hard_sphere_PDF():

    dHS = 1.0

    test_hard_sphere_PDF( V, xHS, dHS, rmax=5.0, N=101 ):
    N = 301

    dHS = 1.0
    V = 1.3
    V = 3.
    lapfun0 = lambda s, V=V, xHS=1.0, dHS=dHS:\
        np.squeeze( hsmix.hard_sphere_LT_PDF( s, V, np.array([xHS]),
                                             np.array([dHS]) ) )

    delt = np.linspace(0.01,6,N)
    ynuminv = hsmix.inv_laplace_euler( lapfun0, delt, tshft=dHS )

    fpack = np.pi/6*dHS**3/V
    lam0 = 2*np.pi/(1-fpack)
    lam1 = np.pi**2*(dHS**2/V)/(1-fpack)**2

    fpack = np.pi/6*dHS**3/V
    zeta = fpack*dHS**3
    ((1-zeta) + 1.5*fpack*dHS**3)/(1.0-zeta)**2

    gij_contact = hsmix.hard_sphere_contact_PDF( V, np.array([xHS]),
                                                np.array([dHS]) )


    r = np.linspace(dHS, 6*dHS,100)
    gij = hsmix.hard_sphere_PDF( r, V, np.array([xHS]), np.array([dHS]) )

    gii =
    gij = 1.0/(2*np.pi)*(lam0 + 0.5*lam1*dHS + 1.0/18*lam1**2/lam0*dHS**2)
    # lapfun = lambda s: np.exp(s*tshft)*lapfun0(s)
    # ynuminv = hsmix.nlinvsteh( lapfun, delt, n=10 )

    plt.plot( delt+dHS, ynuminv/(delt+dHS), 'k-')

