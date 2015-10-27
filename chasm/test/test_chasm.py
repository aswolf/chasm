import numpy as np
import matplotlib.pyplot as plt
import pdb
import hsmix

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

def test_ideal_mix():
    kT = 1.0
    xHS = np.array([.5,.5])
    Fmix, Smix = hsmix.ideal_mix( kT, xHS )

    assert Smix == np.log(2), 'Smix of 50/50 mix should equal log(2)'

    Fmix, Smix = hsmix.ideal_mix( kT, np.array([0.0,1.0]) )
    assert Smix==0, 'Purely 1 component yields Smix=0'


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



