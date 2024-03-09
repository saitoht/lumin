import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import scipy.optimize as optimize
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
import sys, os, pathlib
import cls_subs as subs

class emod:
    """ Class: calculate the elastic moduli using first-principles codes """
    def __init__(self):
        """ Constructor: read parameters from input file """
        subs.get_prms()















### ----------------------------------------------------------------------------- ###
def qjob_dis(ndir: str, alat: float, plat: float, mat: str, ncore: int, exepath: str):
    """ submit job of distorted crystal """
    
    plat = alat * plat
    plat_str1: str = " {plat11} {plat12} {plat13} ".format(plat11=plat[0,0], plat12=plat[0,1], plat13=plat[0,2])
    plat_str2: str = " {plat21} {plat22} {plat23} ".format(plat21=plat[1,0], plat22=plat[1,1], plat23=plat[1,2])
    plat_str3: str = " {plat31} {plat32} {plat33} ".format(plat31=plat[2,0], plat32=plat[2,1], plat33=plat[2,2])
    sub.run(["mkdir "+ndir], shell=True)
    os.chdir(ndir)
    sub.run(["cp ../../"+mat+".scf.in0 "+mat+".scf.in"], shell=True)
    sub.run(["sed -i -e 's/{alat}/"+str(alat)+"/g' "+mat+".scf.in"], shell=True)
    sub.run(["sed -i -e 's/{plat1}/"+str(plat_str1)+"/g' "+mat+".scf.in"], shell=True)
    sub.run(["sed -i -e 's/{plat2}/"+str(plat_str2)+"/g' "+mat+".scf.in"], shell=True)
    sub.run(["sed -i -e 's/{plat3}/"+str(plat_str3)+"/g' "+mat+".scf.in"], shell=True)
    sub.run(["mpirun -np "+str(ncore)+" "+exepath+"/pw.x < "+mat+".scf.in >& "+mat+".scf.out"], shell=True)
    sub.run(["grep ! "+mat+".scf.out > tote.tmp"], shell=True)
    te: float = float(np.loadtxt("tote.tmp",dtype="str",unpack=True,ndmin=0)[4])
    sub.run(["rm work/"+mat+".save/wfc*"], shell=True)
    os.chdir("../")
    volume: float = abs(np.dot(plat[0],np.cross(plat[1],plat[2])))
    return ( te, volume )

### ----------------------------------------------------------------------------- ###
def Bmod(mat: str, nc: int, exe: str, ratio: float = 0.02, ndiv: int = 15):
    """ calculate volume dependencen of the total energy """
    
    print("* --- Bulk Modulus --- *")
    p0 = pathlib.Path("Bmod")
    if ( not p0.is_dir() ):
        sub.run(["mkdir Bmod"], shell=True)
    os.chdir("Bmod")
    delta: float = np.linspace(1.-ratio, 1.+ratio, ndiv)
    (alat0, plat, elements, nelems, natm, pos, volume0) = get_POSCAR("../POSCAR0")
    alat: float = alat0 * delta
    fn: str = "DATA_Bmod"
    p = pathlib.Path(fn)
    if ( not p.exists() ):
        with open(fn,"w") as f:
            f.write("### delta  alat(Bohr)  volume(Bohr^3)  tote(Ry) \n")
    for i, al in enumerate(alat):
        ndir: str = "alat"+str(round(al,6))
        p = pathlib.Path(ndir)
        if ( not p.is_dir() ):
            (te, vol) = qjob_dis(ndir, al, plat, mat, nc, exe)
            string = " {delta}   {alat}   {volume}   {tote} \n".format(delta=delta[i],alat=al,volume=vol,tote=te)
            with open(fn,"a") as f:
                f.write(string)
    os.chdir("../")

### ----------------------------------------------------------------------------- ###
def Emod(mat: str, nc: int, exe: str, epsilon: float, dfmat: float, sym: str):
    """ calculate strain dependence of the total energy """

    print("* --- Elastic Constant "+sym+" --- *")
    p0 = pathlib.Path(sym)
    if ( not p0.is_dir() ):
        sub.run(["mkdir "+sym], shell=True)
    os.chdir(sym)
    fn: str = "DATA_"+sym
    p = pathlib.Path(fn)
    if ( not p.exists() ):
        with open(fn,"w") as f:
            f.write("### delta  alat(Bohr)  volume(Bohr^3)  tote(Ry) \n")
    (alat, plat0, elements, nelems, natm, pos, volume0) = get_POSCAR("../POSCAR0")
    for i, ep in enumerate(epsilon):
        ndir: str = "delta"+str(round(ep,6))
        p = pathlib.Path(ndir)
        if ( not p.is_dir() ):
            plat: float = np.dot(dfmat[i],plat0)
            (tote, volume) = qjob_dis(ndir, alat, plat, mat, nc, exe)
            string = " {delta}   {alat}   {volume}   {tote} \n".format(delta=ep,alat=alat,volume=volume,tote=tote)
            with open(fn,"a") as f:
                f.write(string)
    os.chdir("../")

### ----------------------------------------------------------------------------- ###
def Emod_fit(para: float, sym: str, ndiv: int = 15):
    """ Fitting the elastic modulus """

    def BM_eq(x: float, E0: float, V0: float, B0: float, B0p: float):
        """ Birch-Murnaghan equation of state """
        q: float = (V0/x)**(2./3.) - 1.
        y: float = E0 + (9.*V0*B0/16.) * ((q**3.)*B0p + (q**2.)*(-4.*q+2.))
        return y

    def DeltaE(delta: float, coeff0: float, coeff1: float, coeff2: float):
        return coeff0 + coeff1 * delta + coeff2 * delta**2.

    fn: str = sym+"/DATA_"+sym
    pf = pathlib.Path(fn)
    if ( not pf.exists() ):
        print("*** ERROR: "+fn+" does not exist!!!")
        sys.exit()
    (alat0, plat0, elements, nelems, natm, pos, volume0) = get_POSCAR("POSCAR0")
    delta, alat, volume, tote = np.loadtxt(fn,dtype="float",unpack=True,ndmin=0)
    if ( sym == "Bmod" ):
        cf = optimize.curve_fit(f=BM_eq, xdata=volume, ydata=tote-min(tote), p0=para)
        pm_str: str = "[E0, V0, B0, B0p]"
    else:
        yd: float = np.array([(tote[i]-min(tote))/volume[i] for i in range(len(tote))])
        cf = optimize.curve_fit(f=DeltaE, xdata=delta, ydata=yd, p0=para)
        pm_str: str = "[0th, 1st, 2nd]"
    Emod: float = cf[0]
    print("* Optimized values")
    print(cf)
    print("* parameters "+sym+" = "+pm_str)
    print(Emod)

    if ( sym == "Bmod" ):
        x_cont: float = np.linspace(min(volume)-1.,max(volume)+1,10*ndiv)
        x_plt: float = volume
        fit: float = BM_eq(x_cont, Emod[0], Emod[1], Emod[2], Emod[3])
        fit_comp: float = BM_eq(volume, Emod[0], Emod[1], Emod[2], Emod[3])
        y_plt: float = tote - min(tote)
        plt.xlabel(r"Volume (Bohr$^3$)")
        plt.ylabel(r"$\Delta$Energy (Ry)")
    else:
        x_cont: float = np.linspace(min(delta)-1.e-3,max(delta)+1.e-3,10*ndiv)
        x_plt: float = delta
        fit: float = DeltaE(x_cont, Emod[0], Emod[1], Emod[2])
        fit_comp: float = DeltaE(delta, Emod[0], Emod[1], Emod[2])
        y_plt: float = np.array([(tote[i] - min(tote))/volume[i] for i in range(len(tote))])
        plt.xlabel(r"Delta")
        plt.ylabel(r"$\Delta$Energy/Volume (Ry/Bohr$^3$)")
    r2: float = r2_score(y_plt, fit_comp)
    print("* R2: ", r2)
    print("*")
    plt.plot(x_cont, fit, label="Fit",color="mediumblue")
    plt.scatter(x_plt, y_plt, label="Data",color="white",edgecolor="crimson",s=30)
    plt.legend()
    plt.savefig(sym+"_fit.pdf")
    plt.show()
    return ( Emod )

### ----------------------------------------------------------------------------- ###
def get_Bmod( Bmod0: float, mat: str, nc: int, exe: str, sw_fit: bool ):
    """ get Bulk modulus by fitting against Birch-Murnaghan equation of state """

    if ( abs(Bmod0) < 1.e-5 ):
        Bmod0: float = 1.e2
    (alat, plat, elements, nelems, natm, pos, volume) = get_POSCAR("POSCAR0")
    para: float = [0.0,volume,Bmod0/AU2GPa,3.]
    Bmod(mat, nc, exe)
    Bm: float = 0.0
    if ( sw_fit ):
        Emod: float = Emod_fit(para,"Bmod")
        Bm: float = AU2GPa * Emod[2]
        print("* Bulk modulus by fitting BM_eq (GPa): ", Bm)
    print("* --- Finish calculating Bulk modulus--- *")
    print("*")
    return ( Bm )

### ----------------------------------------------------------------------------- ###
def get_Econst(brav: str, mat: str, nc: int, exe: str, sw_fit: bool, ratio: float = 0.02, ndiv: int = 15):
    """ check symmetry of the system & get elastic constants"""

    para: float = [0.,0.,0.1]
    print("* --- Check Symmetry --- *")
    ep: float = np.linspace(-ratio, ratio, ndiv)
    if ( brav == "cub" ):  # simple cubic
        ### see for example, M. Jamal, S. J. Asadabadi, I. Ahmad, H. A. R. Aliabad,
        ### Elastic constants of cubic crystals, Computational Materials Science 95 (2014) 592-599
        print("*** Cubic system ***")
        print("* There are 3 independent elastic constants: ")
        print("* C11, C12, C44")
        Econst: float = np.array([0.,0.,0.])
        Bmod: float = 0.0
        Gmod: float = 0.0
        Ymod: float = 0.0
        nu: float = 0.0
        # C11 - C12
        sym: str = "cub-uni"
        dfmat: float = np.array([np.diag([1.+d,1.-d,1./(1.-d**2.)]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E10: float = Emod_fit(para, sym)
            E1: float = E10[2]
        # 1.5 * ( C11 + 2*C12 )
        sym: str = "cub-all"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E20: float = Emod_fit(para, sym)
            E2: float = E20[2]
        # 2 * C44
        sym: str = "cub-xy"
        dfmat: float = np.array([[[1.,d,0.],[d,1.,0.],[0.,0.,1./(1.-d**2.)]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E30: float = Emod_fit(para, sym)
            E3: float = E30[2]

        if ( sw_fit ):
            C11: float = AU2GPa * (6.*E1 + 2.*E2)/9.
            C12: float = AU2GPa * (2.*E2 - 3.*E1)/9.
            C44: float = AU2GPa * 0.5*E3
            Econst: float = np.array([C11, C12, C44])
            Bmod: float = (C11 + 2.*C12)/3.
            Cp: float = 0.5*(C11-C12)
            Cpp: float = C12 - C44
            GV: float = 0.2*(3.*C44 + 2.*Cp)
            Ymod: float = 9.*Bmod*GV/(3.*Bmod+GV)
            nu: float = 0.5 - Ymod/(6.*Bmod)
            zeta: float = (C11 + 8.*C12)/(7.*C11 + 2.*C12)
            GR: float = (5.*(C11-C12)*C44)/(4.*C44+3.*(C11-C12))
            Gmod: float = 0.5*(GV + GR)
            Aniso: float = 2.*C44 / (C11-C12)
            lamd: float = Ymod * nu / ((1.+nu)*(1.-2.*nu))
            mu: float = 0.5 * Ymod / (1.+nu)
            print("* C11 (GPa): ", C11)
            print("* C12 (GPa): ", C12)
            print("* C44 (GPa): ", C44)
            print("* Bulk modulus (GPa): ", Bmod)
            print("* C'=0.5(C11-C12) (GPa): ", Cp)
            print("* C''=C12-C44 (GPa): ", Cpp)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Young modulus (GPa): ", Ymod)
            print("* Poisson ratio nu: ", nu)
            print("* Kleinman parameter zeta: ", zeta)
            print("* Reuss shear modulus (GPa): ", GR)
            print("* Hill shear modulus (GPa): ", Gmod)
            print("* Anisotropy constatn: ", Aniso)
            print("* Lames coefficient lambda: ", lamd)
            print("* Lames coefficient mu: ", mu)

    elif ( brav == "tet" ):  # tetragonal
        ### See "A. H. Reshak, M. Jamal, DFT Calculation for Elastic Constants of Tetragonal Strucrure of
        ### Crystalline Solids with WIEN2k Code: A New Package (Tetra-elastic), Int. J. Electrochem. Sci., 8 (2013) 12252."
        print("*** Tetragonal system ***")
        print("* There are 6 independent elastic constants: ")
        print("* C11, C12, C13, C33, C44, C66")
        Econst: float = np.array([0.,0.,0.,0.,0.,0.])
        Bmod: float = 0.0
        Gmod: float = 0.0
        Ymod: float = 0.0
        nu: float = 0.0
        # C11+C12
        sym: str = "tet_pl"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E10: float = Emod_fit(para, sym)
            E1: float = E10[2]
        # C11+C12+2*C33-4*C13
        sym: str = "tet_pl2"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1./((1.+d)**2.)]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E20: float = Emod_fit(para, sym)
            E2: float = E20[2]
        # 0.5*C33
        sym: str = "tet_ax"
        dfmat: float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E30: float = Emod_fit(para, sym)
            E3: float = E30[2]
        # C11-C12
        sym: str = "tet_ortho"
        dfmat: float = np.array([np.diag([np.sqrt((1.+d)/(1.-d)),np.sqrt((1.-d)/(1.+d)),1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E40: float = Emod_fit(para, sym)
            E4: float = E40[2]
        # 4*C44
        sym: str = "tet_yz"
        dfmat: float = np.array([[[1.,0.,d],[0.,1.,d],[d,d,1.+d**2.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E50: float = Emod_fit(para, sym)
            E5: float = E50[2]
        # 2*C66
        sym: str = "tet_xy"
        dfmat: float = np.array([[[np.sqrt(1.+d**2.),d,0.],[d,np.sqrt(1.+d**2.),0.],[0.,0.,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E60: float = Emod_fit(para, sym)
            E6: float = E60[2]

        if ( sw_fit ):
            C11: float = AU2GPa * 0.5 * (E1 + E4)
            C12: float = AU2GPa * 0.5 * (E4 - E1)
            C33: float = AU2GPa * 2. * E3
            C44: float = AU2GPa * 0.25 * E5
            C66: float = AU2GPa * 0.5 * E6
            C13: float = -AU2GPa * 0.25 *(E2 - E1 - 4.*E3)
            Econst: float = np.array([C11,C12,C13,C33,C44,C66])
            Cmod_sq: float = (C11+C12)*C33 - 4.*C13**2.
            BR: float = Cmod_sq/(C11+C12+2.*C33-4.*C13)
            BV: float = (2.*(C11+C12)+4.*C13+C33)/9.
            GR: float = 15./(18.*BV/Cmod_sq+6./(C11-C12)+6./C44+3./C66)
            GV: float = ((C11+C12)+2.*C33-4.*C13+3.*C11-3.*C12+12.*C44+6.*C66)/30.
            Bmod: float = 0.5*(BV+BR)
            Gmod: float = 0.5*(GV+GR)
            Ymod: float = 9.*Bmod*Gmod/(3.*Bmod+Gmod)
            nu: float = 0.5*(1.-Ymod/(3.*Bmod))
            print("* C11 (GPa): ", C11)
            print("* C12 (GPa): ", C12)
            print("* C13 (GPa): ", C13)
            print("* C33 (GPa): ", C33)
            print("* C44 (GPa): ", C44)
            print("* C66 (GPa): ", C66)
            print("* Voigt Bulk modulus (GPa): ", BV)
            print("* Reuss Bulk modulus (GPa): ", BR)
            print("* Hill Bulk modulus (GPa): ", Bmod)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Reuss Shear modulus (GPa):", GR)
            print("* Hill Shear modulus (GPa)", Gmod)
            print("* Young modulus (GPa): ", Ymod)
            print("* Poisson ratio nu: ", nu)
        
    elif ( brav == "hex" ):  # simple hexagonal
        ### see for example, Z. Zhang, Z. H. Fu, R. F. Zhang, D. Legut, and H. B. Guo,
        ### Anomalous mechanical strengths and shear deformation paths of Al2O3 polymorphs with high ionicity, RCS Advances
        print("*** Hexagonal system ***")
        print("* There are 5 independent elastic constants: ")
        print("* C11, C12, C13, C33, C44")
        Econst: float = np.array([0.,0.,0.,0.,0.])
        Bmod: float = 0.0
        Gmod: float = 0.0
        Ymod: float = 0.0
        nu: float = 0.0
        # C11 + C12
        sym: str = "hex-pl"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E10: float = Emod_fit(para, sym)
            E1: float = E10[2]
        # 0.25 * (C11 - C12)
        sym: str = "hex-xy"
        dfmat: float = np.array([[[1.,0.5*d,0.],[0.5*d,1.,0.],[0.,0.,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E20: float = Emod_fit(para, sym)
            E2: float = E20[2]
        # 0.5 * C33
        sym: str = "hex-ax"
        dfmat: float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E30: float = Emod_fit(para, sym)
            E3: float = E30[2]
        # 0.5 * C44
        sym: str = "hex-yz"
        dfmat: float = np.array([[[1.,0.,0.],[0.,1.,0.5*d],[0.,0.5*d,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E40: float = Emod_fit(para, sym)
            E4: float = E40[2]
        # C11 + C12 + 2.*C13 + 0.5*C33
        sym: str = "hex-all"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E50: float = Emod_fit(para, sym)
            E5: float = E50[2]

        if ( sw_fit ):
            C11: float = AU2GPa * 0.5 * (E1 + 4.*E2)
            C12: float = AU2GPa * 0.5 * (4.*E2 - E1)
            C33: float = AU2GPa * 2. * E3
            C44: float = AU2GPa * 2. * E4
            C13: float = AU2GPa * 0.5 * (E5 - E1 - E3)
            C66: float = AU2GPa * 0.5 * (C11-C12)
            Econst: float = np.array([C11,C12,C13,C33,C44])
            BV: float = (2.*(C11+C12) + C33 + 4.*C13)/9.
            GV: float = (C11+C12+2.*C33-4.*C13+12.*C44+12.*C66)
            BR: float = ((C11+C12)*C33 - 2.*C13**2.)/(C11+C12+2.*C33-4.*C13)
            GR: float = (5.*(((C11+C12)*C33-2.*C13**2.)*C44*C66))/(2.*(3.*BV*C44*C66+((C11+C12)*C33-2.*C13**2.)*(C44+C66)))
            Bmod: float = 0.5*(BV+BR)
            Gmod: float = 0.5*(GV+GR)
            Ymod: float = 9.*Bmod*Gmod/(3.*Bmod+Gmod)
            nu: float = 0.5*(1.-Ymod/(3.*Bmod))
            print("* C11 (GPa): ", C11)
            print("* C12 (GPa): ", C12)
            print("* C33 (GPa): ", C33)
            print("* C44 (GPa): ", C44)
            print("* C13 (GPa): ", C13)
            print("* Voigt Bulk modulus (GPa): ", BV)
            print("* Reuss Bulk modulus (GPa): ", BR)
            print("* Hill Bulk modulus (GPa): ", Bmod)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Reuss Shear modulus (GPa): ", GR)
            print("* Hill Shear modulus (GPa): ", Gmod)
            print("* Young modulus (GPa): ", Ymod)
            print("* Possion ratio: ", nu)
        
    elif ( brav == "mono" ):  # Monoclinic
        print("*** Monoclinic system ***")
        print("* There are 8 independent elastic constants: ")
        print("* C11, C22, C33, C44, C55, C66, C12, C13")
        Econst: float = np.array([0.,0.,0.,0.,0.,0.,0.,0.])
        Bmod: float = 0.0
        Gmod: float = 0.0
        Ymod: float = 0.0
        nu: float = 0.0
        # 0.5 * C11
        sym: str = "mono-unix"
        dfmat: float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E10: float = Emod_fit(para, sym)
            E1: float = E10[2]
        # 0.5 * C22
        sym: str = "mono-uniy"
        dfmat: float = np.array([np.diag([1.,1.+d,1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E20: float = Emod_fit(para, sym)
            E2: float = E20[2]
        # 0.5 * C33
        sym: str = "mono-uniz"
        dfmat: float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E30: float = Emod_fit(para, sym)
            E3: float = E30[2]
        # 0.5 * C44
        sym: str = "mono-yz"
        dfmat: float = np.array([[[1.,0.,0.],[0.,1.,0.5*d],[0.,0.5*d,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E40: float = Emod_fit(para, sym)
            E4: float = E40[2]
        # 0.5 * C55
        sym: str = "mono-zx"
        dfmat: float = np.array([[[1.,0.,0.5*d],[0.,1.,0.],[0.5*d,0.,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E50: float = Emod_fit(para, sym)
            E5: float = E50[2]
        # 0.5 * C66
        sym: str = "mono-xy"
        dfmat: float = np.array([[[1.,0.5*d,0.],[0.5*d,1.,0.],[0.,0.,1.]] for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E60: float = Emod_fit(para, sym)
            E6: float = E60[2]
        # 0.5*C11 + C12 + 0.5*C22
        sym: str = "mono-plxy"
        dfmat: float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E70: float = Emod_fit(para, sym)
            E7: float = E70[2]
        # 0.5*C11 + C13 + 0.5*C33
        sym: str = "mono-plzx"
        dfmat: float = np.array([np.diag([1.+d,1.,1.+d]).tolist() for d in ep])
        Emod(mat, nc, exe, ep, dfmat, sym)
        if ( sw_fit ):
            E80: float = Emod_fit(para, sym)
            E8: float = E80[2]

        if ( sw_fit ):
            C11: float = AU2GPa * 2. * E1
            C22: float = AU2GPa * 2. * E2
            C33: float = AU2GPa * 2. * E3
            C44: float = AU2GPa * 2. * E4
            C55: float = AU2GPa * 2. * E5
            C66: float = AU2GPa * 2. * E6
            C12: float = AU2GPa * (E7 - E1 - E2)
            C13: float = AU2GPa * (E8 - E1 - E3)
            Econst: float = np.array([C11,C12,C13,C22,C33,C44,C55,C66])
            print("* C11 (GPa): ", C11)
            print("* C22 (GPa): ", C22)
            print("* C33 (GPa): ", C33)
            print("* C44 (GPa): ", C44)
            print("* C55 (GPa): ", C55)
            print("* C66 (GPa): ", C66)
            print("* C12 (GPa): ", C12)
            print("* C13 (GPa): ", C13)
        
    else:
        print("*** ERROR: brav should be 'cub', or 'hex'!!!")
        print("*** Other crystal systems are not implemented!!!")
        sys.exit()

    print("* --- Finish calculating Elastic constants --- *")
    print("*")

    return ( Econst, Bmod, Gmod, Ymod, nu )

### ----------------------------------------------------------------------------- ###
def Debye_temp(Bmod: float, Gmod: float):
    """ calculate Debye temperature """
    ### See "H. Alipour, A. Hamedani, & G. Alahyarizadeh, First-principles calculations to investigate the thermal response
    ### of the ZrC(1-x)Nx ceramics at extreme conditions, High Temp. Mater. Proc. 42 (2023) 20220241."
    
    print("* --- Debye temperature --- *")
    (alat, plat, elements, nelems, natm, pos, Vol) = get_POSCAR("POSCAR0")
    Mass_ele: float = [ELEMS_MASS[mat] for mat in elements]
    Mass: float = sum([Mass_ele[i] for i, ne in enumerate(nelems) for j in range(ne)])
    rho: float = Mass*uatm/(Vol*((Bohr*1.e-10)**3.))
    n: float = float(natm)
    vl: float = np.sqrt((3.*Bmod + 4.*Gmod)*1.e9/(3.*rho))
    vt: float = np.sqrt(Gmod*1.e9/rho)
    vm: float = ((1./3.)*(2./(vt**3.) + 1./(vl**3.)))**(-1./3.)
    ThetaD: float = (hconst/kB) * (((3.*n*NA*rho)/(4.*np.pi*Mass*1.e-3))**(1./3.)) * vm

    print("* n: ", n)
    print("* Mass: ", Mass)
    print("* rho (kg/m^3): ", rho)
    print("* longitudinal velocity (m/s): ", vl)
    print("* transverse velocity (m/s): ", vt)
    print("* mean velocity (m/s): ", vm)
    print("* Debye temperature (K): ", ThetaD)
    print("* --- Finish --- *")
    print("*")
    
    return ( ThetaD )

### ----------------------------------------------------------------------------- ###
def run_debye():
    """ execute debye temperature program """

    if ( sw_Bmod ):
        p = pathlib.Path(mat+".scf_Bmod.in")
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Bmod.in {mat}.scf.in0".format(mat=mat)], shell=True)
        Bmod0: float = get_Bmod(300., mat, nc, exe, sw_plt)
    if ( sw_Debye ):
        p = pathlib.Path(mat+".scf_Emod.in")
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Emod.in {mat}.scf.in0".format(mat=mat)], shell=True)
        ( Ec, Bmod, Gmod, Ymod, nu ) = get_Econst(brav, mat, nc, exe, sw_plt)
        if ( sw_plt ):
            ThetaD = Debye_temp(Bmod, Gmod)

### ----------------------------------------------------------------------------- ###
if __name__ == "__main__":
    run_debye()
    sys.exit()
