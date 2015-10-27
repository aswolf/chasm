import numpy as np
import matplotlib.pyplot as plt
import pdb
import chasm

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
        iFgas, iSgas = chasm.ideal_gas( V, kT, xHS, mHS )
        Fgas_a[ind] = iFgas
        Sgas_a[ind] = iSgas

    P = -np.diff(Fgas_a)/np.diff(V_a)*chasm.GPA_EV_ANG3

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
        iFgas, iSgas = chasm.ideal_gas( V, kT, xHS, mHS )
        Fgas_a[ind] = iFgas
        Sgas_a[ind] = iSgas

    S = -np.diff(Fgas_a)/np.diff(kT_a)

    assert np.abs( np.log( np.mean(Sgas_a)/S ) ) < TOL

    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

def test_ideal_mix():
    kT = 1.0
    xHS = np.array([.5,.5])
    Fmix, Smix = chasm.ideal_mix( kT, xHS )

    assert Smix == np.log(2), 'Smix of 50/50 mix should equal log(2)'

    Fmix, Smix = chasm.ideal_mix( kT, np.array([0.0,1.0]) )
    assert Smix==0, 'Purely 1 component yields Smix=0'


