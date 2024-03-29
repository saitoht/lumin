import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import scipy.optimize as optimize
from sklearn.metrics import r2_score
import sys, os, pathlib
import cls_subs as subs

prms = subs.subs()

class emod:
    """ Class: calculate the elastic moduli using first-principles codes """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor of emod """
        
        print("* --- Start Elastic Moduli calculation--- *")
        print("*")
        p = pathlib.Path("{mat}.scf_Bmod.in".format(mat=prms.mat))
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Bmod.in {mat}.scf.in0".format(mat=prms.mat)], shell=True)
        if ( prms.sw_Bmod ):
            emod.get_Bmod(self)
        p = pathlib.Path("{mat}.scf_Emod.in".format(mat=prms.mat))
        if ( p.exists() ):
            sub.run(["cp {mat}.scf_Emod.in {mat}.scf.in0".format(mat=prms.mat)], shell=True)
        emod.get_Econst(self)
        emod.calc_Debye_temp(self)
        if ( prms.sw_egap ):
            emod.get_eig(self)
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
                (te, vol) = emod.qjob_dis(self, ndir, al, plat)
                string = " {delta:.10f}   {alat:.10f}   {volume:.10f}   {tote:.10f} \n".format(delta=delta[i],alat=al,volume=vol,tote=te)
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
                plat:float = plat0.dot(dfmat[i])
                (tote, volume) = emod.qjob_dis(self, ndir, alat, plat)
                string = " {delta:.10f}   {alat:.10f}   {volume:.10f}   {tote:.10f} \n".format(delta=ep,alat=alat,volume=volume,tote=tote)
                with open(fn,"a") as f:
                    f.write(string)
        os.chdir("../")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def Emod_fit(self, para:float, sym:str, sw_3rd:bool = False):
        """ Fitting the elastic modulus """

        def BM_eq(x:float, E0:float, V0:float, B0:float, B0p:float):
            """ Birch-Murnaghan equation of state """
            q:float = (V0/x)**(2./3.) - 1.
            y:float = E0 + (9.*V0*B0/16.) * ((q**3.)*B0p + (q**2.)*(-4.*q+2.))
            return y

        def DeltaE(delta:float, coeff0:float, coeff1:float, coeff2:float):
            return coeff0 + coeff1 * delta + coeff2 * delta**2.

        def DeltaE_3rd(delta:float, coeff0:float, coeff1:float, coeff2:float, coeff3:float):
            return coeff0 + coeff1 * delta + coeff2 * delta**2. + coeff3 * delta**3.

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
            if ( not sw_3rd ):
                yd:float = np.array([(tote[i]-min(tote))/volume0 for i in range(len(tote))])
                cf = optimize.curve_fit(f=DeltaE, xdata=delta, ydata=yd, p0=para)
                pm_str:str = "[0th, 1st, 2nd]"
            else:
                yd:float = np.array([(tote[i]-min(tote))/volume0 for i in range(len(tote))])
                cf = optimize.curve_fit(f=DeltaE_3rd, xdata=delta, ydata=yd, p0=para)
                pm_str:str = "[0th, 1st, 2nd, 3rd]"
        Emod:float = cf[0]
        print("* Optimized values")
        print(cf)
        print("* parameters {sym} = {pm_str}".format(sym=sym, pm_str=pm_str))
        print(Emod)

        if ( sym == "Bmod" ):
            x_cont:float = np.linspace(min(volume)-1.,max(volume)+1,10*prms.ndiv_emod)
            x_plt:float = volume
            fit:float = BM_eq(x_cont, Emod[0], Emod[1], Emod[2], Emod[3])
            fit_comp:float = BM_eq(volume, Emod[0], Emod[1], Emod[2], Emod[3])
            y_plt:float = tote - min(tote)
            if ( prms.sw_plt_emod ):
                plt.xlabel(r"Volume (Bohr$^3$)")
                plt.ylabel(r"$\Delta$Energy (Ry)")
        else:
            if ( not sw_3rd ):
                x_cont:float = np.linspace(min(delta)-1.e-3,max(delta)+1.e-3,10*prms.ndiv_emod)
                x_plt:float = delta
                fit:float = DeltaE(x_cont, Emod[0], Emod[1], Emod[2])
                fit_comp:float = DeltaE(delta, Emod[0], Emod[1], Emod[2])
                y_plt:float = np.array([(tote[i] - min(tote))/volume0 for i in range(len(tote))])
            else:
                x_cont:float = np.linspace(min(delta)-1.e-3,max(delta)+1.e-3,10*prms.ndiv_emod)
                x_plt:float = delta
                fit:float = DeltaE_3rd(x_cont, Emod[0], Emod[1], Emod[2], Emod[3])
                fit_comp:float = DeltaE_3rd(delta, Emod[0], Emod[1], Emod[2], Emod[3])
                y_plt:float = np.array([(tote[i] - min(tote))/volume0 for i in range(len(tote))])
            if ( prms.sw_plt_emod ):
                plt.xlabel(r"Delta")
                plt.ylabel(r"$\Delta$Energy/Volume (Ry/Bohr$^3$)")
        r2:float = r2_score(y_plt, fit_comp)
        print("* R2: ", r2)
        print("*")
        if ( prms.sw_plt_emod ):
            plt.plot(x_cont, fit, label="Fit",color="mediumblue")
            plt.scatter(x_plt, y_plt, label="Data",color="white",edgecolor="crimson",s=30)
            plt.legend()
            plt.savefig(sym+"_fit.pdf")
            plt.show()
        return ( Emod )

    ### ----------------------------------------------------------------------------- ###
    def get_Bmod(self):
        """ get Bulk modulus by fitting against Birch-Murnaghan equation of state """

        print("*** Bulk modulus ***")
        Bmod0:float = 1.e2
        (alat, plat, elements, nelems, natm, pos, volume) = prms.get_POSCAR("POSCAR0")
        para:float = [0.0,volume,Bmod0/prms.AU2GPa,3.]
        emod.calc_Bmod(self)
        Emod:float = emod.Emod_fit(self,para,"Bmod")
        self.Bmod:float = prms.AU2GPa * Emod[2]
        print("* Bulk modulus by fitting BM_eq (GPa): ", self.Bmod)
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def get_Econst(self):
        """ check symmetry of the system & get elastic constants"""

        para:float = [0.,0.,0.1]
        para_3rd:float = [0.,0.,0.1,0.]
        ep:float = np.linspace(-prms.dratio, prms.dratio, prms.ndiv_emod)
        if ( prms.brav == "cub" ):  # simple cubic
            """ see for example, M. Jamal, S. J. Asadabadi, I. Ahmad, H. A. R. Aliabad,
             Elastic constants of cubic crystals, Computational Materials Science 95 (2014) 592-599 """
            print("*** Cubic system ***")
            print("* There are 3 independent elastic constants: ")
            print("* C11, C12, C44")
            """ C11 - C12 """
            sym:str = "cub-uni"
            dfmat:float = np.array([np.diag([1.+d,1.-d,1./(1.-d**2.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para, sym)
            E1:float = E10[2]
            """ 2 * C44 """
            sym:str = "cub-xy"
            dfmat:float = np.array([[[1.,d,0.],[d,1.,0.],[0.,0.,1./(1.-d**2.)]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para, sym)
            E2:float = E20[2]
            """ 1.5 * (C11 + 2*C12) """
            sym:str = "cub-vol"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E3:float = E30[2]

            self.C11:float = prms.AU2GPa * 2. * ( 3.*E1 + E3 ) / 9.
            self.C12:float = prms.AU2GPa * ( 2.*E3 - 3.*E1 ) / 9.
            self.C44:float = prms.AU2GPa * 0.5 * E2
            Econst:float = np.array([self.C11, self.C12, self.C44])
            self.Bmod:float = (self.C11 + 2.*self.C12) / 3.
            Cp:float = 0.5*(self.C11-self.C12)
            Cpp:float = self.C12 - self.C44
            GV:float = 0.2*(3.*self.C44 + 2.*Cp)
            self.Ymod:float = 9.*self.Bmod*GV/(3.*self.Bmod+GV)
            self.nu:float = 0.5 - self.Ymod/(6.*self.Bmod)
            zeta:float = (self.C11 + 8.*self.C12)/(7.*self.C11 + 2.*self.C12)
            GR:float = (5.*(self.C11-self.C12)*self.C44)/(4.*self.C44+3.*(self.C11-self.C12))
            self.Gmod:float = 0.5*(GV + GR)
            Aniso:float = 2.*self.C44 / (self.C11-self.C12)
            lamd:float = self.Ymod * self.nu / ((1.+self.nu)*(1.-2.*self.nu))
            mu:float = 0.5 * self.Ymod / (1.+self.nu)
            print("* Calculated values by DFT *")
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

        elif ( prms.brav == "tet" ):  # tetragonal
            """ See A. H. Reshak, M. Jamal, DFT Calculation for Elastic Constants of Tetragonal Strucrure of
                Crystalline Solids with WIEN2k Code: A New Package (Tetra-elastic), Int. J. Electrochem. Sci., 8 (2013) 12252. """
            print("*** Tetragonal system ***")
            print("* There are 6 independent elastic constants: ")
            print("* C11, C12, C13, C33, C44, C66")
            """ C11+C12 """
            sym:str = "tet_pl"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ C11+C12+2*C33-4*C13 """
            sym:str = "tet_pl2"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1./((1.+d)**2.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para, sym)
            E2:float = E20[2]
            """ 0.5*C33 """
            sym:str = "tet_ax"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E3:float = E30[2]
            """ C11-C12 """
            sym:str = "tet_ortho"
            dfmat:float = np.array([np.diag([np.sqrt((1.+d)/(1.-d)),np.sqrt((1.-d)/(1.+d)),1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E40:float = emod.Emod_fit(self, para, sym)
            E4:float = E40[2]
            """ 4*C44 """
            sym:str = "tet_yz"
            dfmat:float = np.array([[[1.,0.,d],[0.,1.,d],[d,d,1.+d**2.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E50:float = emod.Emod_fit(self, para, sym)
            E5:float = E50[2]
            """ 2*C66 """
            sym:str = "tet_xy"
            dfmat:float = np.array([[[np.sqrt(1.+d**2.),d,0.],[d,np.sqrt(1.+d**2.),0.],[0.,0.,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E60:float = emod.Emod_fit(self, para, sym)
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
            print("* Calculated values by DFT *")
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
        
        elif ( prms.brav == "hex" ):  # simple hexagonal
            """ see for example, Z. Zhang, Z. H. Fu, R. F. Zhang, D. Legut, and H. B. Guo,
                Anomalous mechanical strengths and shear deformation paths of Al2O3 polymorphs with high ionicity, RCS Advances """
            print("*** Hexagonal system ***")
            print("* There are 5 independent elastic constants: ")
            print("* C11, C12, C13, C33, C44")
            """ C11 + C12 """
            sym:str = "hex-pp"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ C11 - C12 """
            sym:str = "hex-pm"
            dfmat:float = np.array([np.diag([1.+d,1.-d,1./(1.-d**2.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para, sym)
            E2:float = E20[2]
            """ 0.5*C33 """
            sym:str = "hex-ax"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E3:float = E30[2]
            """ 2.*C44 """
            sym:str = "hex-yz"
            dfmat:float = np.array([[[1./(1.-d**2.),0.,0.],[0.,1.,d],[0.,d,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E40:float = emod.Emod_fit(self, para, sym)
            E4:float = E40[2]
            """ C11 + C12 + 2.*C13 + 0.5*C33 """
            sym:str = "hex-vol"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E50:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E5:float = E50[2]

            self.C11:float = prms.AU2GPa * 0.5 * (E1 + E2)
            self.C12:float = prms.AU2GPa * 0.5 * (E1 - E2)
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 0.5 * E4
            self.C13:float = prms.AU2GPa * 0.5 * (E5 - E1 - E3)
            self.C66:float = 0.5 * (self.C11-self.C12)
            Econst:float = np.array([self.C11,self.C12,self.C13,self.C33,self.C44])
            BV:float = (2.*(self.C11+self.C12) + self.C33 + 4.*self.C13)/9.
            GV:float = (self.C11+self.C12+2.*self.C33-4.*self.C13+12.*self.C44+12.*self.C66)/30.
            BR:float = ((self.C11+self.C12)*self.C33 - 2.*self.C13**2.)/(self.C11+self.C12+2.*self.C33-4.*self.C13)
            GR:float = (5.*(((self.C11+self.C12)*self.C33-2.*self.C13**2.)*self.C44*self.C66))/(2.*(3.*BV*self.C44*self.C66+((self.C11+self.C12)*self.C33-2.*self.C13**2.)*(self.C44+self.C66)))
            self.Bmod:float = 0.5*(BV+BR)
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = 0.5*(1.-self.Ymod/(3.*self.Bmod))
            print("* Calculated values by DFT *")
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
            
        elif ( prms.brav == "mono" ):  # Monoclinic
            print("*** Monoclinic system ***")
            print("* There are 8 independent elastic constants: ")
            print("* C11, C22, C33, C44, C55, C66, C12, C13")
            """ 0.5 * C11 """
            sym:str = "mono-unix"
            dfmat:float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ 0.5 * C22 """
            sym:str = "mono-uniy"
            dfmat:float = np.array([np.diag([1.,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E2:float = E20[2]
            """ 0.5 * C33 """
            sym:str = "mono-uniz"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E3:float = E30[2]
            """ 0.5 * C44 """
            sym:str = "mono-yz"
            dfmat:float = np.array([[[1.,0.,0.],[0.,1.,0.5*d],[0.,0.5*d,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E40:float = emod.Emod_fit(self, para, sym)
            E4:float = E40[2]
            """ 0.5 * C55 """
            sym:str = "mono-zx"
            dfmat:float = np.array([[[1.,0.,0.5*d],[0.,1.,0.],[0.5*d,0.,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E50:float = emod.Emod_fit(self, para, sym)
            E5:float = E50[2]
            """ 0.5 * C66 """
            sym:str = "mono-xy"
            dfmat:float = np.array([[[1.,0.5*d,0.],[0.5*d,1.,0.],[0.,0.,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E60:float = emod.Emod_fit(self, para, sym)
            E6:float = E60[2]
            """ 0.5*C11 + C12 + 0.5*C22 """
            sym:str = "mono-plxy"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E70:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E7:float = E70[2]
            """ 0.5*C11 + C13 + 0.5*C33 """
            sym:str = "mono-plzx"
            dfmat:float = np.array([np.diag([1.+d,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E80:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
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
            print("* Calculated values by DFT *")
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
        
        """ See Y. Zhou, A. M. Tehrani, A. O. Oliynyk, A. C. Duke, & J. Brogch, Identifying an efficient,
                thermally robust inorganic phosphor host via machine learning, Nature Commn. 9 (2018) 4377. """
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

    ### ----------------------------------------------------------------------------- ###
    def get_eig(self, path:str = ".", code: str = "ecalj"):
        """ get eigenvalues from bnd***.spin* in ecalj and EIGENVAL in VASP """
    
        self.eigu:float = [ [] for i in range(prms.nkpath*prms.nkpt) ]
        self.eigd:float = [ [] for i in range(prms.nkpath*prms.nkpt) ]
        self.ehomo:float = -1000.0
        self.elumo:float = 1000.0
    
        if ( code == "ecalj" ):
            for k in range(prms.nkpath):
                for isp in range(prms.nspin):
                    if ( k+1 < 10 ):
                        fn:str = path+"/bnd00"+str(k+1)+".spin"+str(isp+1)
                    elif ( k+1 > 9 and k+1 < 100 ):
                        fn:str = path+"/bnd0"+str(k+1)+".spin"+str(isp+1)
                    else:
                        fn:str = path+"/bnd"+str(k+1)+".spin"+str(isp+1)
                    p = pathlib.Path(fn)
                    if ( not p.exists() ):
                        print("*** ERROR in m_dat.get_eig: "+fn+" doesn't exist!!!")
                        sys.exit()
                    bnd_org:float = np.loadtxt(fname=fn, comments="#", unpack=False, ndmin=0)
                    count = 1
                    count_bnd = 0
                    for bnd in bnd_org:
                        nb = int(bnd[0])
                        ikp = ( count + count_bnd ) % ( prms.nkpt+1 ) + k * prms.nkpt - 1
                        if ( nb > prms.nbndmin-1 and nb < prms.nbndmax+1 ):
                            bd = bnd[2]
                            if ( bd <= 0.0 and bd >= self.ehomo ):
                                self.ehomo = bd
                            elif ( bd >= 0.0 and bd <= self.elumo ):
                                self.elumo = bd

                            if ( isp == 0 ):
                                self.eigu[ikp].append(bnd[2])
                            else:
                                self.eigd[ikp].append(bnd[2])
                            count += 1
                        else:
                            count += 1
                        if ( count % prms.nkpt == 0 ):
                            count_bnd += 1

        elif ( code == "VASP" ):
            fn:str = path + "/EIGENVAL" 
            p = pathlib.Path(fn)
            if ( not p.exists() ):
                print("*** ERROR in m_dat.get_eig: "+fn+" doesn't exist!!!")
                sys.exit()

            with open(file=fn, mode="r") as f:
                data:str = f.readlines()
            
            count:int = 0
            count_kp:int = 0
            sw_skip:bool = False
            for dat in data:
                if ( count < 6 ):
                    count += 1
                elif ( ( count-count_kp-6 ) % ( prms.nbnds+1 ) == 0 ):
                    if ( sw_skip ):
                        sw_skip = False
                    else:
                        sw_skip = True
                        count_kp += 1
                    count += 1
                else:
                    nb = ( count-count_kp-6 ) % ( prms.nbnds+1 )
                    if ( nb > prms.nbndmin-1 and nb < prms.nbndmax+1 ):
                        if ( prms.nspin == 1 ):
                            nn, eu, occu = dat.split()
                        else:
                            nn, eu, ed, occu, occd = dat.split()
                            print(eu, ed)
                            self.eigd[count_kp-1].append(ed)
                        self.eigu[count_kp-1].append(eu)
                    count += 1

        elif ( code == "QE" ):
            for isp in range(prms.nspin):
                if ( isp == 0 ):
                    fn:str = path+"/"+prms.mat+".bands_up.dat"
                else:
                    fn:str = path+"/"+prms.mat+".bands_dn.dat"
                p = pathlib.Path(fn)
                if ( not p.exists() ):
                    print("*** ERROR in m_dat.get_eig: "+fn+" doesn't exist!!!")
                    sys.exit()

            sub.run(["grep Fermi "+prms.mat+".scf.out > grep.out"], shell=True)
            with open("grep.out", "r") as f:
                ef:float = float(f.readlines()[0].split()[4])
            with open(file=fn, mode="r") as f:
                data:str = f.readlines()

            for i, line in enumerate(data):
                if ( i == 0 ):
                    nbnd:int = int(line.split()[2].replace(",",""))
                    nkpt:int = int(line.split()[4])
                    ncyc:int = prms.nbnd // 10
                    nres:int = prms.nbnd % 10
                    if ( not nres == 0 ):
                        ncyc += 1
                    ikp:int = -1
                else:
                    pass
                if ( (i-1)//(ncyc+1) ):
                    ikp += 1
                    nb:int = 0
                else:
                    pass
                if ( isp == 0 ):
                    for bnd in line.split().replace("\n",""):
                        if ( nb > prms.nbndmin-1 and nb < prms.nbndmax+1 ):
                            self.eigu[ikp].append(float(bnd)-prms.ef)
                            nb += 1
                            if ( float(bnd) <= ef and float(bnd) >= self.ehomo ):
                                self.ehomo = float(bnd)
                            elif ( float(bnd) >= ef and float(bnd) <= self.elumo ):
                                self.elumo = float(bnd)
                            else:
                                pass
                        else:
                            pass
                else:
                    for bnd in line.split().replace("\n",""):
                        if ( nb > prms.nbndmin-1 and nb < prms.nbndmax+1 ):
                            self.eigd[ikp].append(float(bnd)-self.ef)
                            nb += 1
                            if ( float(bnd) <= self.ef and float(bnd) >= self.ehomo ):
                                self.ehomo = float(bnd)
                            elif ( float(bnd) >= self.ef and float(bnd) <= self.elumo ):
                                self.elumo = float(bnd)
                            else:
                                pass
                        else:
                            pass
            
        else:
            print("*** ERROR in m_dat.get_eig: code should be 'ecalj', 'QE' or 'VASP'!!!")
            sys.exit()

        # E_g = LUMO - HOMO
        self.egap: float = self.elumo - self.ehomo
        print("* Energy gap (eV): ", self.egap)
