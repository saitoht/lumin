import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import scipy.optimize as optimize
from sklearn.metrics import r2_score
import sys, os, pathlib
import cls_subs as subs

class emod:
    """ Class: calculate the elastic moduli using first-principles codes """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor of reading parameters from input file """
        
        prms = subs.subs()
        print("* --- Start Elastic Moduli calculation--- *")
        print("*")
        p = pathlib.Path("{mat}.scf_Bmod.in".format(mat=prms.mat))
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Bmod.in {mat}.scf.in0".format(mat=prms.mat)], shell=True)
        emod.get_Bmod(self)
        p = pathlib.Path("{mat}.scf_Emod.in".format(mat=prms.mat))
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Emod.in {mat}.scf.in0".format(mat=prms.mat)], shell=True)
        emod.get_Econst(self)
        emod.calc_Debye_temp(self)
        print("* --- Finish Elastic Moduli calculation--- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def qjob_dis(self, ndir:str, alat:float, plat:float):
        """ execute a job of distorted crystal """
        
        plat = alat * plat
        plat_str1:str = " {plat11} {plat12} {plat13} ".format(plat11=plat[0,0], plat12=plat[0,1], plat13=plat[0,2])
        plat_str2:str = " {plat21} {plat22} {plat23} ".format(plat21=plat[1,0], plat22=plat[1,1], plat23=plat[1,2])
        plat_str3:str = " {plat31} {plat32} {plat33} ".format(plat31=plat[2,0], plat32=plat[2,1], plat33=plat[2,2])
        sub.run(["mkdir "+ndir], shell=True)
        os.chdir(ndir)
        sub.run(["cp ../../{mat}.scf.in0 {mat}.scf.in".format(mat=prms.mat)], shell=True)
        sub.run(["sed -i -e 's/{alat}/"+str(alat)+"/g' "+prms.mat+".scf.in"], shell=True)
        sub.run(["sed -i -e 's/{plat1}/"+str(plat_str1)+"/g' "+prms.mat+".scf.in"], shell=True)
        sub.run(["sed -i -e 's/{plat2}/"+str(plat_str2)+"/g' "+prms.mat+".scf.in"], shell=True)
        sub.run(["sed -i -e 's/{plat3}/"+str(plat_str3)+"/g' "+prms.mat+".scf.in"], shell=True)
        sub.run(["mpirun -np {nc} {exe}/pw.x < {mat}.scf.in >& {mat}.scf.out".format(nc=prms.nc,exe=prms.exe,mat=prms.mat)], shell=True)
        sub.run(["grep ! {mat}.scf.out > tote.tmp".format(mat=prms.mat)], shell=True)
        te:float = float(np.loadtxt("tote.tmp",dtype="str",unpack=True,ndmin=0)[4])
        sub.run(["rm work/{mat}.save/wfc*".format(mat=prms.mat)], shell=True)
        os.chdir("../")
        volume:float = abs(np.dot(plat[0],np.cross(plat[1],plat[2])))
        return ( te, volume )

    ### ----------------------------------------------------------------------------- ###
    def calc_Bmod(self):
        """ calculate volume dependencen of the total energy """
        
        print("*** First-principles for Bulk Modulus ***")
        p = pathlib.Path("Bmod")
        if ( not p.is_dir() ):
            sub.run(["mkdir Bmod"], shell=True)
        os.chdir("Bmod")
        delta:float = np.linspace(1.-prms.dratio, 1.+prms.dratio, prms.ndiv_emod)
        (alat0, plat, elements, nelems, natm, pos, volume0) = prms.get_POSCAR("../POSCAR0")
        alat:float = alat0 * delta
        fn:str = "DATA_Bmod"
        p = pathlib.Path(fn)
        if ( not p.exists() ):
            with open(fn,"w") as f:
                f.write("### delta  alat(Bohr)  volume(Bohr^3)  tote(Ry) \n")
            for i, al in enumerate(alat):
                ndir:str = "alat"+str(round(al,6))
                p = pathlib.Path(ndir)
            if ( not p.is_dir() ):
                (te, vol) = qjob_dis(ndir, al, plat)
                string = " {delta}   {alat}   {volume}   {tote} \n".format(delta=delta[i],alat=al,volume=vol,tote=te)
                with open(fn,"a") as f:
                    f.write(string)
        os.chdir("../")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def calc_Emod(self, epsilon:float, dfmat:float, sym:str):
        """ calculate strain dependence of the total energy """
        
        print("*** First-principles for Elastic Moduli {sym} ***".format(sym=sym))
        p = pathlib.Path(sym)
        if ( not p.is_dir() ):
            sub.run(["mkdir "+sym], shell=True)
        os.chdir(sym)
        fn:str = "DATA_"+sym
        p = pathlib.Path(fn)
        if ( not p.exists() ):
            with open(fn,"w") as f:
                f.write("### delta  alat(Bohr)  volume(Bohr^3)  tote(Ry) \n")
        (alat, plat0, elements, nelems, natm, pos, volume0) = prms.get_POSCAR("../POSCAR0")
        for i, ep in enumerate(epsilon):
            ndir:str = "delta"+str(round(ep,6))
            p = pathlib.Path(ndir)
            if ( not p.is_dir() ):
                plat:float = np.dot(dfmat[i],plat0)
                (tote, volume) = qjob_dis(ndir, alat, plat)
                string = " {delta}   {alat}   {volume}   {tote} \n".format(delta=ep,alat=alat,volume=volume,tote=tote)
                with open(fn,"a") as f:
                    f.write(string)
        os.chdir("../")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def Emod_fit(self, para:float, sym:str, ndiv:int = 15):
        """ Fitting the elastic modulus """
        
        def BM_eq(x:float, E0:float, V0:float, B0:float, B0p:float):
            """ Birch-Murnaghan equation of state """
            q:float = (V0/x)**(2./3.) - 1.
            y:float = E0 + (9.*V0*B0/16.) * ((q**3.)*B0p + (q**2.)*(-4.*q+2.))
            return y

        def DeltaE(delta:float, coeff0:float, coeff1:float, coeff2:float):
            return coeff0 + coeff1 * delta + coeff2 * delta**2.

        fn:str = "{sym}/DATA_{sym}".format(sym=sym)
        pf = pathlib.Path(fn)
        if ( not pf.exists() ):
            print("*** ERROR in emod.Emod_fit: {fn} does not exist!!!".format(fn=fn))
            sys.exit()
        (alat0, plat0, elements, nelems, natm, pos, volume0) = prms.get_POSCAR("POSCAR0")
        delta, alat, volume, tote = np.loadtxt(fn,dtype="float",unpack=True,ndmin=0)
        if ( sym == "Bmod" ):
            cf = optimize.curve_fit(f=BM_eq, xdata=volume, ydata=tote-min(tote), p0=para)
            pm_str:str = "[E0, V0, B0, B0p]"
        else:
            yd:float = np.array([(tote[i]-min(tote))/volume[i] for i in range(len(tote))])
            cf = optimize.curve_fit(f=DeltaE, xdata=delta, ydata=yd, p0=para)
            pm_str:str = "[0th, 1st, 2nd]"
        Emod:float = cf[0]
        print("* Optimized values")
        print(cf)
        print("* parameters {sym} = {pm_str}".format(sym=sym, pm_str=pm_str))
        print(Emod)

        if ( sym == "Bmod" ):
            x_cont:float = np.linspace(min(volume)-1.,max(volume)+1,10*ndiv)
            x_plt:float = volume
            fit:float = BM_eq(x_cont, Emod[0], Emod[1], Emod[2], Emod[3])
            fit_comp:float = BM_eq(volume, Emod[0], Emod[1], Emod[2], Emod[3])
            y_plt:float = tote - min(tote)
            plt.xlabel(r"Volume (Bohr$^3$)")
            plt.ylabel(r"$\Delta$Energy (Ry)")
        else:
            x_cont:float = np.linspace(min(delta)-1.e-3,max(delta)+1.e-3,10*ndiv)
            x_plt:float = delta
            fit:float = DeltaE(x_cont, Emod[0], Emod[1], Emod[2])
            fit_comp:float = DeltaE(delta, Emod[0], Emod[1], Emod[2])
            y_plt:float = np.array([(tote[i] - min(tote))/volume[i] for i in range(len(tote))])
            plt.xlabel(r"Delta")
            plt.ylabel(r"$\Delta$Energy/Volume (Ry/Bohr$^3$)")
        r2:float = r2_score(y_plt, fit_comp)
        print("* R2: ", r2)
        print("*")
        plt.plot(x_cont, fit, label="Fit",color="mediumblue")
        plt.scatter(x_plt, y_plt, label="Data",color="white",edgecolor="crimson",s=30)
        plt.legend()
        plt.savefig(sym+"_fit.pdf")
        if ( prms.sw_plt_emod ):
            plt.show()
        return ( Emod )

    ### ----------------------------------------------------------------------------- ###
    def get_Bmod(self, Bmod0:float):
        """ get Bulk modulus by fitting against Birch-Murnaghan equation of state """

        print("*** Bulk modulus ***")
        if ( abs(Bmod0) < 1.e-5 ):
            Bmod0:float = 1.e2
        (alat, plat, elements, nelems, natm, pos, volume) = prms.get_POSCAR("POSCAR0")
        para:float = [0.0,volume,Bmod0/prms.AU2GPa,3.]
        ccd.calc_Bmod(self)
        Emod:float = ccd.Emod_fit(self,para,"Bmod")
        self.Bmod:float = prms.AU2GPa * Emod[2]
        print("* Bulk modulus by fitting BM_eq (GPa): ", self.Bmod)
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def get_Econst(self):
        """ check symmetry of the system & get elastic constants"""
        
        para:float = [0.,0.,0.1]
        ep:float = np.linspace(-prms.dratio, prms.dratio, prms.ndiv_emod)
        if ( prms.brav == "cub" ):  # simple cubic
            ### see for example, M. Jamal, S. J. Asadabadi, I. Ahmad, H. A. R. Aliabad,
            ### Elastic constants of cubic crystals, Computational Materials Science 95 (2014) 592-599
            print("*** Cubic system ***")
            print("* There are 3 independent elastic constants: ")
            print("* C11, C12, C44")
            # C11 - C12
            sym:str = "cub-uni"
            dfmat:float = np.array([np.diag([1.+d,1.-d,1./(1.-d**2.)]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E10:float = ccd.Emod_fit(self, para, sym)
            E1:float = E10[2]
            # 2 * C44
            sym:str = "cub-xy"
            dfmat:float = np.array([[[1.,d,0.],[d,1.,0.],[0.,0.,1./(1.-d**2.)]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E20:float = ccd.Emod_fit(self, para, sym)
            E2:float = E20[2]

            self.C11:float = prms.AU2GPa * (6.*E1 + 9.*self.Bmod)/9.
            self.C12:float = prms.AU2GPa * (9.*self.Bmod - 3.*E1)/9.
            self.C44:float = prms.AU2GPa * 0.5*E2
            Econst:float = np.array([self.C11, self.C12, self.C44])
            Cp:float = 0.5*(self.C11-self.C12)
            Cpp:float = self.C12 - self.C44
            GV:float = 0.2*(3.*self.C44 + 2.*self.Cp)
            self.Ymod:float = 9.*self.Bmod*GV/(3.*self.Bmod+GV)
            self.nu:float = 0.5 - self.Ymod/(6.*self.Bmod)
            zeta:float = (self.C11 + 8.*self.C12)/(7.*self.C11 + 2.*self.C12)
            GR:float = (5.*(self.C11-self.C12)*self.C44)/(4.*self.C44+3.*(self.C11-self.C12))
            self.Gmod:float = 0.5*(GV + GR)
            Aniso:float = 2.*self.C44 / (self.C11-self.C12)
            lamd:float = self.Ymod * self.nu / ((1.+self.nu)*(1.-2.*self.nu))
            mu:float = 0.5 * self.Ymod / (1.+self.nu)
            print("* C11 (GPa): ", self.C11)
            print("* C12 (GPa): ", self.C12)
            print("* C44 (GPa): ", self.C44)
            print("* Bulk modulus (GPa): ", self.Bmod)
            print("* C'=0.5(C11-C12) (GPa): ", Cp)
            print("* C''=C12-C44 (GPa): ", Cpp)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Young modulus (GPa): ", self.Ymod)
            print("* Poisson ratio nu: ", self.nu)
            print("* Kleinman parameter zeta: ", zeta)
            print("* Reuss shear modulus (GPa): ", GR)
            print("* Hill shear modulus (GPa): ", self.Gmod)
            print("* Anisotropy constatn: ", Aniso)
            print("* Lames coefficient lambda: ", lamd)
            print("* Lames coefficient mu: ", mu)
            print("*")

        elif ( brav == "tet" ):  # tetragonal
            ### See "A. H. Reshak, M. Jamal, DFT Calculation for Elastic Constants of Tetragonal Strucrure of
            ### Crystalline Solids with WIEN2k Code: A New Package (Tetra-elastic), Int. J. Electrochem. Sci., 8 (2013) 12252."
            print("*** Tetragonal system ***")
            print("* There are 6 independent elastic constants: ")
            print("* C11, C12, C13, C33, C44, C66")
            # C11+C12
            sym:str = "tet_pl"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E10:float = ccd.Emod_fit(self, para, sym)
            E1:float = E10[2]
            # C11+C12+2*C33-4*C13
            sym:str = "tet_pl2"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1./((1.+d)**2.)]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E20:float = ccd.Emod_fit(self, para, sym)
            E2:float = E20[2]
            # 0.5*C33
            sym:str = "tet_ax"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E30:float = ccd.Emod_fit(self, para, sym)
            E3:float = E30[2]
            # C11-C12
            sym:str = "tet_ortho"
            dfmat:float = np.array([np.diag([np.sqrt((1.+d)/(1.-d)),np.sqrt((1.-d)/(1.+d)),1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E40:float = ccd.Emod_fit(self, para, sym)
            E4:float = E40[2]
            # 4*C44
            sym:str = "tet_yz"
            dfmat:float = np.array([[[1.,0.,d],[0.,1.,d],[d,d,1.+d**2.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E50:float = ccd.Emod_fit(self, para, sym)
            E5:float = E50[2]
            # 2*C66
            sym:str = "tet_xy"
            dfmat:float = np.array([[[np.sqrt(1.+d**2.),d,0.],[d,np.sqrt(1.+d**2.),0.],[0.,0.,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E60:float = ccd.Emod_fit(self, para, sym)
            E6:float = E60[2]

            self.C11:float = prms.AU2GPa * 0.5 * (E1 + E4)
            self.C12:float = prms.AU2GPa * 0.5 * (E4 - E1)
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 0.25 * E5
            self.C66:float = prms.AU2GPa * 0.5 * E6
            self.C13:float = -prms.AU2GPa * 0.25 *(E2 - E1 - 4.*E3)
            Econst:float = np.array([self.C11,self.C12,self.C13,self.C33,self.C44,self.C66])
            Cmod_sq:float = (self.C11+self.C12)*self.C33 - 4.*self.C13**2.
            BR:float = Cmod_sq/(self.C11+self.C12+2.*self.C33-4.*self.C13)
            BV:float = (2.*(self.C11+self.C12)+4.*self.C13+self.C33)/9.
            GR:float = 15./(18.*BV/Cmod_sq+6./(self.C11-self.C12)+6./self.C44+3./self.C66)
            GV:float = ((self.C11+self.C12)+2.*self.C33-4.*self.C13+3.*self.C11-3.*self.C12+12.*self.C44+6.*self.C66)/30.
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = 0.5*(1.-self.Ymod/(3.*self.Bmod))
            print("* C11 (GPa): ", self.C11)
            print("* C12 (GPa): ", self.C12)
            print("* C13 (GPa): ", self.C13)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C66 (GPa): ", self.C66)
            print("* Voigt Bulk modulus (GPa): ", BV)
            print("* Reuss Bulk modulus (GPa): ", BR)
            print("* Hill Bulk modulus (GPa): ", self.Bmod)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Reuss Shear modulus (GPa):", GR)
            print("* Hill Shear modulus (GPa)", self.Gmod)
            print("* Young modulus (GPa): ", self.Ymod)
            print("* Poisson ratio nu: ", self.nu)
            print("*")
        
        elif ( brav == "hex" ):  # simple hexagonal
            ### see for example, Z. Zhang, Z. H. Fu, R. F. Zhang, D. Legut, and H. B. Guo,
            ### Anomalous mechanical strengths and shear deformation paths of Al2O3 polymorphs with high ionicity, RCS Advances
            print("*** Hexagonal system ***")
            print("* There are 5 independent elastic constants: ")
            print("* C11, C12, C13, C33, C44")
            # C11 + C12
            sym:str = "hex-pl"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E10:float = ccd.Emod_fit(self, para, sym)
            E1:float = E10[2]
            # 0.25 * (C11 - C12)
            sym:str = "hex-xy"
            dfmat:float = np.array([[[1.,0.5*d,0.],[0.5*d,1.,0.],[0.,0.,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E20:float = ccd.Emod_fit(self, para, sym)
            E2:float = E20[2]
            # 0.5 * C33
            sym:str = "hex-ax"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E30:float = ccd.Emod_fit(self, para, sym)
            E3:float = E30[2]
            # 0.5 * C44
            sym:str = "hex-yz"
            dfmat:float = np.array([[[1.,0.,0.],[0.,1.,0.5*d],[0.,0.5*d,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E40:float = ccd.Emod_fit(self, para, sym)
            E4:float = E40[2]
            # C11 + C12 + 2.*C13 + 0.5*C33
            sym:str = "hex-all"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E50:float = ccd.Emod_fit(self, para, sym)
            E5:float = E50[2]

            self.C11:float = prms.AU2GPa * 0.5 * (E1 + 4.*E2)
            self.C12:float = prms.AU2GPa * 0.5 * (4.*E2 - E1)
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 2. * E4
            self.C13:float = prms.AU2GPa * 0.5 * (E5 - E1 - E3)
            self.C66:float = prms.AU2GPa * 0.5 * (self.C11-self.C12)
            Econst:float = np.array([self.C11,self.C12,self.C13,self.C33,self.C44])
            BV:float = (2.*(self.C11+self.C12) + self.C33 + 4.*self.C13)/9.
            GV:float = (self.C11+self.C12+2.*self.C33-4.*self.C13+12.*self.C44+12.*self.C66)
            BR:float = ((self.C11+self.C12)*self.C33 - 2.*self.C13**2.)/(self.C11+self.C12+2.*self.C33-4.*self.C13)
            GR:float = (5.*(((self.C11+self.C12)*self.C33-2.*self.C13**2.)*self.C44*self.C66))/(2.*(3.*BV*self.C44*self.C66+((self.C11+self.C12)*self.C33-2.*self.C13**2.)*(self.C44+self.C66)))
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = 0.5*(1.-self.Ymod/(3.*self.Bmod))
            print("* C11 (GPa): ", self.C11)
            print("* C12 (GPa): ", self.C12)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C13 (GPa): ", self.C13)
            print("* Voigt Bulk modulus (GPa): ", BV)
            print("* Reuss Bulk modulus (GPa): ", BR)
            print("* Hill Bulk modulus (GPa): ", self.Bmod)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Reuss Shear modulus (GPa): ", GR)
            print("* Hill Shear modulus (GPa): ", self.Gmod)
            print("* Young modulus (GPa): ", self.Ymod)
            print("* Possion ratio: ", self.nu)
            print("*")
            
        elif ( brav == "mono" ):  # Monoclinic
            print("*** Monoclinic system ***")
            print("* There are 8 independent elastic constants: ")
            print("* C11, C22, C33, C44, C55, C66, C12, C13")
            # 0.5 * C11
            sym:str = "mono-unix"
            dfmat:float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E10:float = ccd.Emod_fit(self, para, sym)
            E1:float = E10[2]
            # 0.5 * C22
            sym:str = "mono-uniy"
            dfmat:float = np.array([np.diag([1.,1.+d,1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E20:float = ccd.Emod_fit(self, para, sym)
            E2:float = E20[2]
            # 0.5 * C33
            sym:str = "mono-uniz"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E30:float = ccd.Emod_fit(self, para, sym)
            E3:float = E30[2]
            # 0.5 * C44
            sym:str = "mono-yz"
            dfmat:float = np.array([[[1.,0.,0.],[0.,1.,0.5*d],[0.,0.5*d,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E40:float = ccd.Emod_fit(self, para, sym)
            E4:float = E40[2]
            # 0.5 * C55
            sym:str = "mono-zx"
            dfmat:float = np.array([[[1.,0.,0.5*d],[0.,1.,0.],[0.5*d,0.,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E50:float = ccd.Emod_fit(self, para, sym)
            E5:float = E50[2]
            # 0.5 * C66
            sym:str = "mono-xy"
            dfmat:float = np.array([[[1.,0.5*d,0.],[0.5*d,1.,0.],[0.,0.,1.]] for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E60:float = ccd.Emod_fit(self, para, sym)
            E6:float = E60[2]
            # 0.5*C11 + C12 + 0.5*C22
            sym:str = "mono-plxy"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E70:float = ccd.Emod_fit(self, para, sym)
            E7:float = E70[2]
            # 0.5*C11 + C13 + 0.5*C33
            sym:str = "mono-plzx"
            dfmat:float = np.array([np.diag([1.+d,1.,1.+d]).tolist() for d in ep])
            ccd.calc_Emod(self, ep, dfmat, sym)
            E80:float = ccd.Emod_fit(self, para, sym)
            E8:float = E80[2]

            self.C11:float = prms.AU2GPa * 2. * E1
            self.C22:float = prms.AU2GPa * 2. * E2
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 2. * E4
            self.C55:float = prms.AU2GPa * 2. * E5
            self.C66:float = prms.AU2GPa * 2. * E6
            self.C12:float = prms.AU2GPa * (E7 - E1 - E2)
            self.C13:float = prms.AU2GPa * (E8 - E1 - E3)
            Econst:float = np.array([self.C11,self.C12,self.C13,self.C22,self.C33,self.C44,self.C55,self.C66])
            print("* C11 (GPa): ", self.C11)
            print("* C22 (GPa): ", self.C22)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C55 (GPa): ", self.C55)
            print("* C66 (GPa): ", self.C66)
            print("* C12 (GPa): ", self.C12)
            print("* C13 (GPa): ", self.C13)
            print("*")
            
        else:
            print("*** ERROR in emod.get_Econst: brav should be 'cub', or 'hex'!!!")
            print("*** Other crystal systems are not implemented!!!")
            sys.exit()

    ### ----------------------------------------------------------------------------- ###
    def calc_Debye_temp(self):
        """ calculate Debye temperature """
        
        ### See "H. Alipour, A. Hamedani, & G. Alahyarizadeh, First-principles calculations to investigate the thermal response
        ### of the ZrC(1-x)Nx ceramics at extreme conditions, High Temp. Mater. Proc. 42 (2023) 20220241."
        (alat, plat, elements, nelems, natm, pos, Vol) = prms.get_POSCAR("POSCAR0")
        Mass_ele:float = [prms.ELEMS_MASS[ele] for ele in elements]
        Mass:float = sum([Mass_ele[i] for i, ne in enumerate(nelems) for j in range(ne)])
        self.rho:float = Mass*prms.uatm/(Vol*((prms.Bohr*1.e-10)**3.))
        n:float = float(natm)
        self.vl:float = np.sqrt((3.*self.Bmod + 4.*self.Gmod)*1.e9/(3.*self.rho))
        self.vt:float = np.sqrt(self.Gmod*1.e9/self.rho)
        self.vm:float = ((1./3.)*(2./(self.vt**3.) + 1./(self.vl**3.)))**(-1./3.)
        self.ThetaD:float = (prms.hconst/prms.kB) * (((3.*n*prms.NA*self.rho)/(4.*np.pi*Mass*1.e-3))**(1./3.)) * self.vm

        print("* calc_Debye_temp *")
        print("* n: ", n)
        print("* Mass: ", Mass)
        print("* rho (kg/m^3): ", self.rho)
        print("* longitudinal velocity (m/s): ", self.vl)
        print("* transverse velocity (m/s): ", self.vt)
        print("* mean velocity (m/s): ", self.vm)
        print("* Debye temperature (K): ", self.ThetaD)
        print("*")
