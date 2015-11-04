import numpy as np
import matplotlib.pyplot as plt
import bondpoly
import hsmix


MPROTON_AU = 1836.15
EV_HA = 27.211
ANG_BOHR = .529177
GPA_EV_ANG3 = 160.2176487 # GPa/(eV/Ang^3)
KT300 = 0.025852 # thermal energy in eV at 300K

class chasm( object ):
    """
    components are defined on per-cation basis
    model_d keys:
        comp_name: list of strings
        comp_mass: array of doubles
        comp_noxy: array of doubles

        catoxypot: array of dicts
        oxyoxypot: dict

        generic ionionpot: {'type':'buck','A':A,'rho':rho,'C':C}

    """
    def __init__( self, model_d ):
        """
        oxide components are defined on a per-cation-basis
          (use comp_ncat~=1 to convert to per-cation-basis)
        """

        self.model_d = model_d

        # comp_name, comp_mass, comp_noxy, comp_ncat=1
        # store component info
        self.comp_name = comp_name
        self.comp_mass = 1.0*comp_mass/comp_ncat
        self.comp_noxy = 1.0*comp_noxy/comp_ncat

        # initialize and store bondpoly geometry
        poly_geom = bondpoly.fit_poly_geom()
        self.poly_edgelen_f = poly_geom.edgelen_f
        self.poly_bondang_f = bondpoly.calc_bond_ang

        # mHS_a, noxyHS_a, bondlendev_a, catoxypot_d, oxyoxypot_d
        pass


#====================================================================
# SECT 1: module functions
#====================================================================
def eval_popwt_energy( V, T, X_a, param_d ):
    Fpopwt
    popCN
    return Fpopwt, popCN

def eval_approx_energy( V, T, X_a, param_d ):
    """
    environ conditions: V, T, X_a
    condition-dependent parameters: cn_a, distdev_a
    """
    cn_a = param_d['cn_a']
    distdev_a = param_d['distdev_a']
    distHSadj = self.model_d['distHSadj']

    # determine ideal closest packing geom
    distoxy = self.calc_closest_packing_dist( V, X_a )
    edgelenfac_a = self.edgelen_f( cn_a )
    bondlen_ideal_a = distoxy/edgelenfac_a

    # calc cat-oxy bondlen (based on deviation from ideal values)
    bondlen_a = bondlen_ideal_a * np.exp( distdev_a )

    # calc effective hardsphere size (based on adjustment to cat-oxy bondlen)
    dHS_a = bondlen_a * np.exp( distHSadj )

    Fgas = self.eval_ideal_gas_energy( V, T, X_a )
    Fmix = self.eval_ideal_mix_energy( T, X_a )
    Fhsmix = self.eval_hs_mix_energy( V, T, X_a, dHS_a )
    Fstruc = self.eval_struc_energy( V, T, X_a, bondlen_a )

    Fvib = self.eval_vib_energy( V, T, X_a, param_d )

    Fapprox = Fgas + Fmix + Fhsmix + Fstruc + Fvib
    return Fapprox

def eval_ideal_gas_energy( self, V, T, X_a ):
    kT = KT300*T/300.0
    mass = self.comp_mass
    Fgas, Sgas = hsmix.ideal_gas( V, kT, X_a, mass )
    return Fgas

def eval_ideal_mix_energy( self, T, X_a ):
    kT = KT300*T/300.0
    Fmix = hsmix.ideal_mix( kT, X_a )
    return Fmix

def eval_hs_mix_energy( self, V, T, X_a, dHS ):
    kT = KT300*T/300.0
    Fhsmix_kT, = hsmix.hard_sphere_mix( V, X_a, dHS )
    Fhsmix = kT*Fhsmix_kT
    return Fhsmix


def eval_struc_energy( V, T, X_a, bondlen_a ):
    kT = KT300*T/300.0

    # eval struc energy according to bondlen
    Fstruc = 0.0
    return Fstruc


#====================================================================
# SECT 2: library functions
#====================================================================
def calc_closest_packing_dist( self, V, X_a ):
    noxyHS = self.comp_noxy
    Noxy = np.sum( x_a*noxyHS )
    Voxy = V/Noxy
    fpack_hcp = np.pi/np.sqrt(18)
    distoxy = (6/np.pi*fpack_hcp*Voxy)**(1.0/3)
    return distoxy


def eval_buck_energy( param, cn, dist ):
    A   = param[0]
    rho = param[1]
    C   = param[2]

    Epair = cn*(A*np.exp(-dist/rho) - C/dist**6)
    return Epair
