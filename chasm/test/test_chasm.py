import numpy as np
import matplotlib.pyplot as plt
import pdb
import chasm

def test_ideal_gas():
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
        print V
        iFgas, iSgas = chasm.ideal_gas( V, kT, xHS, mHS )
        Fgas_a[ind] = iFgas
        Sgas_a[ind] = iSgas

    P = -np.diff(Fgas_a)/np.diff(V_a)*chasm.GPA_EV_ANG3

    assert np.abs(np.log(P*1e4/1.013)) < .03, \
        'Must recover 1 bar atmospheric pressure'





