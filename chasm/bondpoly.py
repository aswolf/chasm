import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def fit_poly_geom():
    """
    obtain fitted geometry of bonding polyhedra

    - use weighted fit to emphasize ideal structures
    """
    poly_geom = get_poly_geom()
    cn = poly_geom['cn']
    edgelen = poly_geom['edgelen']
    bondang = poly_geom['bondang']

    wt = np.ones(edgelen.shape)
    wt[cn==2]  = 100
    wt[cn==4]  = 3
    wt[cn==6]  = 3
    wt[cn==12] = 3

    def resid_f( param, cn=cn, edgelen=edgelen, wt=wt):
        resid = wt*(edgelen - decay_model(param,cn))
        return resid

    # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
    param0 = np.array([1,10,1.3,15])
    fitoutput = optimize.leastsq( resid_f, param0, full_output=True)
    paramf = fitoutput[0]
    # covf = fitoutput[1]
    # print fitoutput
    # print '=========='
    # print paramf
    # print np.sqrt(np.diag(covf))

    cnmod = np.linspace(cn[0],cn[-1]+4,10001)
    edgelen_mod = decay_model( paramf, cnmod )
    bondang_mod = calc_bond_ang( edgelen_mod )

    # plt.clf()
    # plt.plot(cn,edgelen,'ko',cnmod,edgelen_mod,'r-')
    # plt.plot(cn,bondang,'ko',cnmod,bondang_mod,'r-')

    def edgelen_model( cn, param=paramf ):
        return decay_model(param, cn)

    poly_geom['edgelen_param'] = paramf
    poly_geom['edgelen_f'] = lambda cn, param=paramf: decay_model(param,cn)

    return poly_geom


def calc_closest_packing_nndist( V, nobj ):
    fpack = np.pi/np.sqrt(18)
    nndist = (6/np.pi*fpack*V/nobj)**(1.0/3)
    return nndist


def calc_bond_ang( edgelen, bondlen=1.0 ):
    ratio = edgelen/bondlen
    bondang = 2*np.arcsin(ratio/2.0)*180/np.pi
    return bondang



def decay_model( param, cn ):
    A0 = param[0]
    cn0 = param[1]
    A1 = param[2]
    cn1 = param[3]

    ymod = A0*np.exp(-cn/cn0) + A1*np.exp(-cn/cn1)
    return ymod


def get_poly_geom():
    poly_geom = {}

    poly_pos = get_poly_pos()

    cn = np.array( poly_pos.keys() )
    bondang = []
    edgelen = []
    for key in poly_pos:
        # print np.sum(poly_pos[key]**2,axis=1)
        ipos = np.expand_dims(poly_pos[key],axis=0)
        jpos = np.expand_dims(poly_pos[key],axis=1)
        distij = np.sqrt(np.sum((ipos-jpos)**2,axis=2))
        np.fill_diagonal(distij,np.nan)
        idistNN = np.mean(np.nanmin(distij,axis=0))
        ibondang = 2*np.arcsin(idistNN/2.0)*180/np.pi
        edgelen.append(idistNN)
        bondang.append(ibondang)

    bondang = np.array( bondang )
    edgelen = np.array( edgelen )

    poly_geom['cn'] = cn
    poly_geom['pos'] = poly_pos
    poly_geom['bondang'] = bondang
    poly_geom['edgelen'] = edgelen

    return poly_geom

def get_poly_pos():
    poly_pos = {}

    poly_pos[2] = np.array([[0, 0,  1.0],
                            [0, 0, -1.0]])
    poly_pos[3] = np.array([[      0, 1.0000, 0.0],
                            [ 0.8660,-0.5000, 0.0],
                            [-0.8660,-0.5000, 0.0]])
    poly_pos[4] = np.array([[+0.720811059899,-0.605940036634,-0.336553246800],
                            [-0.843553324940,-0.170505595776,-0.509259884334],
                            [-0.156171632120,-0.175543065384,+0.972005685949],
                            [+0.278913892764,+0.951988698026,-0.126192548332]])
    poly_pos[5] = np.array([[+0.855565360836,-0.156802100197,+0.493377152604],
                            [-0.658633419068,-0.744679382224,-0.107956643971],
                            [-0.196931925292,+0.901481477843,-0.385420500112],
                            [+0.443792542701,-0.268572623162,-0.854936795986],
                            [-0.443792549436,+0.268572625140,+0.854936791869]])
    poly_pos[6] = np.array([[-0.477973841915,+0.243304922042,+0.844004574250],
                            [-0.169605180501,+0.917224510901,-0.360462590776],
                            [+0.477973836299,-0.243304924988,-0.844004576582],
                            [+0.861844008545,+0.315439234995,+0.397143543271],
                            [-0.861844004834,-0.315439241394,-0.397143546240],
                            [+0.169605180866,-0.917224508446,+0.360462596851]])
    poly_pos[7] = np.array([[-0.570584671282,-0.374416093280,-0.730921146219],
                            [-0.500743916118,-0.546814264721,+0.671006475654],
                            [+0.442735410739,-0.827918915986,-0.344290029816],
                            [+0.079015190293,+0.743083893065,+0.664517063416],
                            [+0.844210175839,-0.137265964975,+0.518138238185],
                            [-0.795376086465,+0.596517091729,-0.107444126667],
                            [+0.500743902593,+0.546814263304,-0.671006486901]])
    poly_pos[8] = np.array([[-0.358740453187,-0.111796186592,-0.926718349830],
                            [+0.230506662846,-0.878717440372,+0.417998012398],
                            [-0.687412312891,+0.726243033921,-0.005947080580],
                            [+0.704082308404,-0.438460053911,-0.558588295722],
                            [-0.832316101440,-0.552053522700,+0.049867979287],
                            [-0.270751592548,+0.183947828588,+0.944910986014],
                            [+0.398985382895,+0.806565799153,-0.436190641660],
                            [+0.815646098192,+0.264270542800,+0.514667390372]])
    poly_pos[9] = np.array([[+0.914109572223,-0.182781178690,-0.361931942064],
                            [+0.293329304506,+0.734642489361,-0.611766566546],
                            [-0.480176899428,-0.046026929940,+0.875963279468],
                            [-0.705684904851,+0.704780196051,-0.072757750931],
                            [+0.370605109670,+0.769162968265,+0.520615194684],
                            [-0.904030464226,-0.412626217894,-0.111662545460],
                            [-0.162180419233,-0.247163999394,-0.955304908927],
                            [+0.063327560246,-0.997971078243,-0.006583851785],
                            [+0.610701141906,-0.322016246902,+0.723429092590]])
    poly_pos[10] = np.array([[+0.978696890330,+0.074682616274,+0.191245663177],
                             [+0.537258145625,+0.448413180814,-0.714338368164],
                             [-0.227939324473,-0.303819959434,-0.925060590777],
                             [+0.274577116268,+0.833436432027,+0.479573895237],
                             [-0.599426405232,+0.240685139624,+0.763386303437],
                             [-0.424664555168,+0.830194107787,-0.361161679833],
                             [-0.402701180119,-0.893328907767,+0.199487398294],
                             [+0.552788606831,-0.770301636525,-0.317899583084],
                             [+0.290107593166,-0.385278374104,+0.876012647646],
                             [-0.978696887344,-0.074682599351,-0.191245685067]])
    poly_pos[11] = np.array([[+0.153486836562,-0.831354332797,+0.534127105044],
                             [+0.092812115769,+0.691598091278,-0.716294626049],
                             [+0.686120068086,+0.724987503180,+0.060269166267],
                             [+0.101393837471,+0.257848797505,+0.960850293931],
                             [-0.143059218646,-0.243142754178,-0.959382958495],
                             [-0.909929380017,+0.200934944687,-0.362841110384],
                             [-0.405338453688,+0.872713317547,+0.272162090194],
                             [+0.896918545883,-0.184616420020,+0.401813264476],
                             [+0.731466092268,-0.415052523977,-0.541007170195],
                             [-0.439821168531,-0.864743799130,-0.242436592901],
                             [-0.773718984882,-0.203685975092,+0.599892453681]])
    poly_pos[12] = np.array([[-0.519657468039,+0.854102387855,-0.021569120792],
                             [+0.519657456743,-0.854102395012,+0.021569109539],
                             [-0.118190661273,+0.429081803333,-0.895499734024],
                             [-0.768920201670,+0.071819589916,+0.635298095360],
                             [+0.521506258421,+0.836678341438,-0.167333724625],
                             [+0.119333282457,+0.615878166459,+0.778751341429],
                             [+0.915718081273,+0.043626880273,+0.399445980012],
                             [+0.768920197428,-0.071819593191,-0.635298100124],
                             [+0.118190681568,-0.429081793598,+0.895499736009],
                             [-0.521506259512,-0.836678336747,+0.167333744677],
                             [-0.915718079511,-0.043626889027,-0.399445983095],
                             [-0.119333287884,-0.615878161696,-0.778751344364]])
    poly_pos[13] = np.array([[-0.754267365796,+0.173499508397,-0.633228759202],
                             [+0.330321246216,-0.372227441281,+0.867372242037],
                             [-0.533035489131,-0.458931921778,+0.710812674690],
                             [-0.686807855664,-0.708384153523,-0.162747843106],
                             [-0.321357401597,+0.941504078299,-0.101486407881],
                             [+0.617869468935,+0.786068318474,-0.018273424674],
                             [+0.019352413289,+0.576222568464,+0.817063666854],
                             [+0.311368579326,-0.948900371701,+0.051358469549],
                             [-0.882505894394,+0.313694805173,+0.350398224264],
                             [-0.038448968078,-0.536051111410,-0.843309482225],
                             [+0.184148321180,+0.468947287620,-0.863815858410],
                             [+0.936881686583,+0.003852303897,+0.349625321023],
                             [+0.813567705110,-0.236010693892,-0.531419365069]])
    poly_pos[14] = np.array([[+0.183504500440,-0.774327848934,-0.605592668948],
                             [+0.153461674965,+0.930458382726,+0.332711154506],
                             [-0.708967620726,-0.280759148895,-0.646946066588],
                             [-0.976939230344,-0.107848809998,+0.184277981311],
                             [+0.372001336525,+0.124628732488,+0.919827529846],
                             [-0.280635908687,-0.474838246429,+0.834129562169],
                             [-0.374884405704,-0.924316474887,+0.071419441413],
                             [+0.609952952660,-0.713378465035,+0.345034144926],
                             [-0.545724577900,+0.497759969260,+0.674106592519],
                             [-0.009781260146,+0.151939133540,-0.988341452459],
                             [+0.836141855120,-0.174860987699,-0.519894636535],
                             [-0.609952964510,+0.713378454079,-0.345034146631],
                             [+0.930390223148,+0.274617378657,+0.242815419632],
                             [+0.421433390542,+0.757547917835,-0.498512837869]])
    poly_pos[15] = np.array([[+0.216001983700,-0.898615365314,+0.381881615504],
                             [-0.444561363553,-0.418726663881,+0.791854263732],
                             [-0.635186145240,-0.770065881729,-0.059473512527],
                             [+0.500548727571,-0.192538217024,+0.844026069688],
                             [-0.221229771469,+0.531169628316,-0.817872981685],
                             [-0.891941312897,+0.079142510044,-0.445170930601],
                             [-0.169315707509,+0.513429353508,+0.841262438331],
                             [-0.263965919566,-0.440254059158,-0.858194824444],
                             [+0.598524035730,+0.006919955250,-0.801074960833],
                             [-0.899631198936,+0.223429384229,+0.375157321885],
                             [+0.985305483557,-0.166404497880,+0.038505157549],
                             [+0.627905121857,+0.627079786181,+0.460983838881],
                             [-0.350119335865,+0.936470382954,+0.020968369105],
                             [+0.527048413266,+0.773124394001,-0.352843650184],
                             [+0.420616988873,-0.804160710609,-0.420008214424]])
    poly_pos[16] = np.array([[+0.733571208539,-0.509367037098,-0.449909439244],
                             [-0.527936513082,+0.197513112343,+0.825997341768],
                             [+0.210490445990,+0.934016014999,-0.288631002964],
                             [-0.394388830898,+0.876999734693,+0.274461136430],
                             [+0.901500079093,+0.427890594987,-0.064863287901],
                             [-0.126116414468,-0.264945505828,-0.955980401966],
                             [-0.807197676003,-0.422051151370,-0.412679945579],
                             [+0.264549795600,-0.885039709623,+0.383038011218],
                             [-0.616921323093,-0.639391147814,+0.458897636963],
                             [-0.969209372877,+0.238904864368,+0.059646100528],
                             [-0.080651424152,-0.925976042468,-0.368868156063],
                             [+0.393234069013,+0.625539858854,+0.673844827799],
                             [+0.874906967098,-0.244744309460,+0.417897142739],
                             [+0.175293851258,-0.219349658857,+0.959769656152],
                             [-0.509883948542,+0.523218324223,-0.682833028065],
                             [+0.478759086918,+0.286782058057,-0.829786591763]])
    return poly_pos
