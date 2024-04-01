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
        
        print("* --- Start Elastic Moduli calculation --- *")
        print("*")
        if ( prms.sw_relax ):
            exejob:str = "relax"
        else:
            exejob:str = "scf"
        p = pathlib.Path("{mat}.{job}_Bmod.in".format(mat=prms.mat,job=exejob))
        if ( p.exists() and prms.code=="qe" ):
            sub.run(["cp {mat}.{job}_Bmod.in {mat}.{job}.in0".format(mat=prms.mat,job=exejob)], shell=True)
        p = pathlib.Path("ctrl_Bmod.{mat}".format(mat=prms.mat))
        if ( p.exists() and prms.code=="ecalj" ):
            sub.run(["cp ctrl_Bmod.{mat} ctrl.{mat}".format(mat=prms.mat)], shell=True)
        if ( prms.sw_Bmod ):
            emod.get_Bmod(self)
        p = pathlib.Path("{mat}.{job}_Emod.in".format(mat=prms.mat,job=exejob))
        if ( p.exists() and prms.code=="qe" ):
            sub.run(["cp {mat}.{job}_Emod.in {mat}.{job}.in0".format(mat=prms.mat,job=exejob)], shell=True)
        p = pathlib.Path("ctrl_Emod.{mat}".format(mat=prms.mat))
        if ( p.exists() and prms.code=="ecalj" ):
            sub.run(["cp ctrl_Emod.{mat} ctrl.{mat}".format(mat=prms.mat)], shell=True)
        emod.get_Econst(self)
        emod.calc_Debye_temp(self)
        if ( prms.sw_egap ):
            emod.get_eig(self)
        print("* --- Finish Elastic Moduli calculation --- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def qjob_dis(self, ndir:str, alat:float, plat:float):
        """ execute a job of distorted crystal """

        skip:bool = False
        plat = alat * plat
        plat_str1:str = " {plat11} {plat12} {plat13} ".format(plat11=plat[0,0], plat12=plat[0,1], plat13=plat[0,2])
        plat_str2:str = " {plat21} {plat22} {plat23} ".format(plat21=plat[1,0], plat22=plat[1,1], plat23=plat[1,2])
        plat_str3:str = " {plat31} {plat32} {plat33} ".format(plat31=plat[2,0], plat32=plat[2,1], plat33=plat[2,2])
        sub.run(["mkdir "+ndir], shell=True)
        os.chdir(ndir)
        if ( prms.code == "qe" ):
            if ( prms.relax ):
                exejob:str = "relax"
            else:
                exejob:str = "scf"
            sub.run(["cp ../../{mat}.{job}.in0 {mat}.{job}.in".format(mat=prms.mat,job=exejob)], shell=True)
            sub.run(["sed -i -e 's/{alat}/"+str(alat)+"/g' "+prms.mat+"."+exejob+".in"], shell=True)
            sub.run(["sed -i -e 's/{plat1}/"+str(plat_str1)+"/g' "+prms.mat+"."+exejob+".in"], shell=True)
            sub.run(["sed -i -e 's/{plat2}/"+str(plat_str2)+"/g' "+prms.mat+"."+exejob+".in"], shell=True)
            sub.run(["sed -i -e 's/{plat3}/"+str(plat_str3)+"/g' "+prms.mat+"."+exejob+".in"], shell=True)
            sub.run(["mpirun -np {nc} {exe}/pw.x < {mat}.{job}.in >& {mat}.{job}.out".format(nc=prms.nc,exe=prms.exe,mat=prms.mat,job=exejob)], shell=True)
            sub.run(["grep ! {mat}.{job}.out | tail -1 > tote.tmp".format(mat=prms.mat,job=exejob)], shell=True)
            te:float = float(np.loadtxt("tote.tmp",dtype="str",unpack=True,ndmin=0)[4])
            sub.run(["rm work/{mat}.save/wfc*".format(mat=prms.mat)], shell=True)
        elif ( prms.code == "ecalj" ):
            sub.run(["cp ../../*.{mat} .".format(mat=prms.mat)], shell=True)
            sub.run(["sed -i -e 's/{alat}/"+str(alat)+"/g' ctrl."+prms.mat], shell=True)
            sub.run(["sed -i -e 's/{plat1}/"+str(plat_str1)+"/g' ctrl."+prms.mat], shell=True)
            sub.run(["sed -i -e 's/{plat2}/"+str(plat_str2)+"/g' ctrl."+prms.mat], shell=True)
            sub.run(["sed -i -e 's/{plat3}/"+str(plat_str3)+"/g' ctrl."+prms.mat], shell=True)
            sub.run(["{exe}/lmfa {mat} >& llmfa".format(exe=prms.exe,mat=prms.mat)], shell=True)
            sub.run(["mpirun -np {nc} {exe}/lmf {mat} >& llmf".format(nc=prms.nc,exe=prms.exe,mat=prms.mat)], shell=True)
            p = pathlib.Path("./save."+prms.mat)
            if ( not p.exists() ):
                print("*** ERROR in cls_dat.qjob: The save file does not exist! Something went wrong in the scf calculation in ecalj!!!")
                skip:bool = True
            else:
                sub.run(["grep c save.{mat} | tail -1 > tote.tmp".format(mat=prms.mat)], shell=True)
                p = pathlib.Path("tote.tmp")
                if ( p.stat().st_size == 0 ):
                    print("*** ERROR: The calculation does not be converged!!!")
                    skip:bool = True
                else:
                    data:str = np.loadtxt("tote.tmp",dtype=str,unpack=True,ndmin=0)
                    if ( prms.nspin == 1 ):
                        if ( data[1][0:8] == "ehk(eV)=" ):
                            te:float = float(data[1].replace("ehk(eV)=",""))/prms.Ry
                        else:
                            te:float = float(data[1].replace("ehf(eV)=",""))/prms.Ry
                    elif ( prms.nspin == 2 ):
                        if ( data[1] == "mmom=" ):
                            te:float = float(data[3].replace("ehk(eV)=",""))/prms.Ry
                        else:
                            te:float = float(data[2].replace("ehk(eV)=",""))/prms.Ry
                    else:
                        pass
        os.chdir("../")
        volume:float = abs(np.dot(plat[0],np.cross(plat[1],plat[2])))
        return ( te, volume, skip )

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
                f.write("### 1:delta  2:alat(Bohr)  3:volume(Bohr^3)  4:tote(Ry) \n")
        for i, al in enumerate(alat):
            ndir:str = "alat"+str(round(al,6))
            p = pathlib.Path(ndir)
            if ( not p.is_dir() ):
                (te, vol, skip) = emod.qjob_dis(self, ndir, al, plat)
                if ( skip ):
                    pass
                else:
                    string = " {delta:.10f}   {alat:.10f}   {volume:.10f}   {tote:.10f} \n".format(delta=delta[i],alat=al,volume=vol,tote=te)
                    with open(fn,"a") as f:
                        f.write(string)
        os.chdir("../")

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
                f.write("### 1:delta  2:alat(Bohr)  3:volume(Bohr^3)  4:tote(Ry) \n")
        (alat, plat0, elements, nelems, natm, pos, volume0) = prms.get_POSCAR("../POSCAR0")
        for i, ep in enumerate(epsilon):
            ndir:str = "delta"+str(round(ep,6))
            p = pathlib.Path(ndir)
            if ( not p.is_dir() ):
                plat:float = plat0.dot(dfmat[i])
                (tote, volume, skip) = emod.qjob_dis(self, ndir, alat, plat)
                if ( skip ):
                    pass
                else:
                    string = " {delta:.10f}   {alat:.10f}   {volume:.10f}   {tote:.10f} \n".format(delta=ep,alat=alat,volume=volume,tote=tote)
                    with open(fn,"a") as f:
                        f.write(string)
        os.chdir("../")

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
        if ( prms.brav == "cub" ):  # cubic
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

        elif ( prms.brav == "ortho" ):  # orthorombic
            """ see P.W.O. Nyawere et al., Physica B 434 (2014) 122-128.,
            and L. Liu, et al., Crystals 7 (2017) 111 """
            print("*** Orthorombic system ***")
            print("* There are 9 independent elastic constatns: ")
            print("* C11, C22, C33, C44, C55, C66, C12, C13, C23")
            """ 0.5*C11 """
            sym:str = "ortho-x"
            dfmat:float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ 0.5*C22 """
            sym:str = "ortho-y"
            dfmat:float = np.array([np.diag([1.,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E2:float = E20[2]
            """ 0.5*C33 """
            sym:str = "ortho-z"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E3:float = E30[2]
            """ 2.*C44 """
            sym:str = "ortho-yz"
            dfmat:float = np.array([[[(1.-d**2.)**(-1./3.),0.,0.],[0.,(1.-d**2.)**(-1./3.),d*(1.-d**2.)**(-1./3.)],[0.,d*(1.-d**2.)**(-1./3.),(1.-d**2.)**(-1./3.)]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E40:float = emod.Emod_fit(self, para, sym)
            E4:float = E40[2]
            """ 2.*C55 """
            sym:str = "ortho-zx"
            dfmat:float = np.array([[[(1.-d**2.)**(-1./3.),0.,d*(1.-d**2.)**(-1./3.)],[0.,(1.-d**2.)**(-1./3.),0.],[d*(1.-d**2.)**(-1./3.),0.,(1.-d**2.)**(-1./3.)]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E50:float = emod.Emod_fit(self, para, sym)
            E5:float = E50[2]
            """ 2.*C66 """
            sym:str = "ortho-xy"
            dfmat:float = np.array([[[(1.-d**2.)**(-1./3.),d*(1.-d**2.)**(-1./3.),0.],[d*(1.-d**2.)**(-1./3.),(1.-d**2.)**(-1./3.),0.],[0.,0.,(1.-d**2.)**(-1./3.)]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E60:float = emod.Emod_fit(self, para, sym)
            E6:float = E60[2]
            """ C11 + C22 - 2.*C12 """
            sym:str = "ortho-pmxy"
            dfmat:float = np.array([np.diag([(1.+d)*(1.-d**2.)**(-1./3.),(1.-d)*(1.-d**2.)**(-1./3.),(1.-d**2.)**(-1./3.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E70:float = emod.Emod_fit(self, para, sym)
            E7:float = E70[2]
            """ C11 + C33 - 2.*C13 """
            sym:str = "ortho-pmzx"
            dfmat:float = np.array([np.diag([(1.+d)*(1.-d**2.)**(-1./3.),(1.-d**2.)**(-1./3.),(1.-d)*(1.-d**2.)**(-1./3.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E80:float = emod.Emod_fit(self, para, sym)
            E8:float = E80[2]
            """ C22 + C33 - 2.*C23 """
            sym:str = "ortho-pmyz"
            dfmat:float = np.array([np.diag([(1.-d**2.)**(-1./3.),(1.+d)*(1.-d**2.)**(-1./3.),(1.+d)*(1.-d**2.)**(-1./3.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E90:float = emod.Emod_fit(self, para, sym)
            E9:float = E90[2]

            self.C11:float = prms.AU2GPa * 2. * E1
            self.C22:float = prms.AU2GPa * 2. * E2
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 0.5 * E4
            self.C55:float = prms.AU2GPa * 0.5 * E5
            self.C66:float = prms.AU2GPa * 0.5 * E6
            self.C12:float = -prms.AU2GPa * 0.5 * (E7 - 2.*E1 - 2.*E2)
            self.C13:float = -prms.AU2GPa * 0.5 * (E8 - 2.*E3 - 2.*E1)
            self.C23:float = -prms.AU2GPa * 0.5 * (E9 - 2.*E2 - 2.*E3)
            Econst:float = np.array([self.C11,self.C22,self.C33,self.C44,self.C55,self.C66,self.C12,self.C13,self.C23])
            BV:float = (self.C11+2.*self.C12+2.*self.C13+self.C22+2.*self.C23+self.C33)/9.
            GV:float = (self.C11-self.C12-self.C13+self.C22-self.C23+self.C33+3.*self.C44+3.*self.C55+3.*self.C66)/15.
            chi:float = self.C13*(self.C12*self.C23-self.C13*self.C22) + self.C23*(self.C12*self.C13-self.C23*self.C11) + self.C33*(self.C11*self.C22-self.C12**2.)
            BR:float = chi/(self.C11*(self.C22+self.C33-2.*self.C23) + self.C22*(self.C33-2.*self.C13) - 2.*self.C33*self.C12 + self.C12*(2.*self.C23-self.C12) + self.C13*(2.*self.C12-self.C13) + self.C23*(2.*self.C13-self.C23))
            GR:float = 15./(4.*(self.C11*(self.C22+self.C33+self.C23) + self.C22*(self.C33+self.C13) + self.C33*self.C12 - self.C12*(self.C23+self.C12) - self.C13*(self.C12+self.C13) - self.C23*(self.C13+self.C23))/chi + 3.*(1./self.C44)+1./self.C55+1./self.C66)
            self.Bmod:float = 0.5*(BV+BR)
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = (3.*self.Bmod-2.*self.Gmod)/(6.*self.Bmod+2.*self.Gmod)
            print("* Calculated values by DFT *")
            print("* C11 (GPa): ", self.C11)
            print("* C22 (GPa): ", self.C22)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C55 (GPa): ", self.C55)
            print("* C66 (GPa): ", self.C66)
            print("* C12 (GPa): ", self.C12)
            print("* C13 (GPa): ", self.C13)
            print("* C23 (GPa): ", self.C23)
            print("* Voigt Bulk modulus (GPa): ", BV)
            print("* Reuss Bulk modulus (GPa): ", BR)
            print("* Hill Bulk modulus (GPa): ", self.Bmod)
            print("* Voigt Shear modulus (GPa): ", GV)
            print("* Reuss Shear modulus (GPa): ", GR)
            print("* Hill Shear modulus (GPa): ", self.Gmod)
            print("* Young modulus (GPa): ", self.Ymod)
            print("* Possion ratio: ", self.nu)
            print("*")
        
        elif ( prms.brav == "hex" ):  # hexagonal
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

        elif ( prms.brav == "trig1" ):  # Trigonal (D3, C3v, D3d)
            print("*** Trigonal (D3, C3v, D3d) system ***")
            print("* There are 6 independent elastic constants: ")
            print("* C11, C33, C12, C13, C44, C24")
            """ 0.5 * C11 """
            sym:str = "trig1-x"
            dfmat:float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ 0.5 * C33 """
            sym:str = "trig1-z"
            dfmat:float = np.array([np.diag([1.,1.,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E2:float = E20[2]
            """ C11 - C12 """
            sym:str = "trig1-pm"
            dfmat:float = np.array([np.diag([1.+d,1.-d,1./(1.-d**2.)]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E30:float = emod.Emod_fit(self, para, sym)
            E3:float = E30[2]
            """ 0.5 * ( 2.*C11 + C33 + 2.*C12 + 4.*C13 ) """
            sym:str = "trig1-vol"
            dfmat:float = np.array([np.diag([1.+d,1.+d,1.+d]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E40:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E4:float = E40[2]
            """ C44 """
            sym:str = "trig1-yz"
            dfmat:float = np.array([[[(1.-d**2.)**(-1.),0.,0.],[0.,1.,d],[0.,d,1.]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E50:float = emod.Emod_fit(self, para, sym)
            E5:float = E50[2]
            """ 0.5 * ( 2.*C11 - 2.*C12 + C44 - 4.*C24 ) """
            sym:str = "trig1-pmyz"
            dfmat:float = np.array([[[1.+d,0.,0.],[0.,1.-d,d],[0.,d,(1.-d**2.)**(-1.)]] for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E60:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E6:float = E60[2]

            self.C11:float = prms.AU2GPa * 2. * E1
            self.C33:float = prms.AU2GPa * 2. * E2
            self.C12:float = -prms.AU2GPa * (E3-2.*E1)
            self.C13:float = prms.AU2GPa * 0.5 * (E4-E1+0.5*E2) - 0.5 * self.C12
            self.C44:float = prms.AU2GPa * E5
            self.C24:float = -prms.AU2GPa * 0.5 * (E6+E1+0.25*E5) - 0.5 * self.C12
            Econst:float = np.array([self.C11,self.C33,self.C12,self.C13,self.C44,self.C24])
            BV:float = (2.*(self.C11+self.C12+2.*self.C13)+self.C33)/9.
            GV:float = ((7.*self.C11-5.*self.C12-4.*self.C13)+(2.*self.C33+12.*self.C44))/30.
            Stiffness:float = np.array([[self.C11,self.C12,self.C13,-self.C24,0.,0.],
                                        [self.C12,self.C11,self.C13,self.C24,0.,0.],
                                        [self.C13,self.C13,self.C33,0.,0.,0.],
                                        [-self.C24,self.C24,0.,self.C44,0.,0.],
                                        [0.,0.,0.,0.,self.C44,-self.C24],
                                        [0.,0.,0.,0.,-self.C24,0.5*(self.C11-self.C12)]])
            Compliance:float = np.linalg.inv(Stiffness)
            S11:float = Compliance[0,0]
            S22:float = Compliance[1,1]
            S33:float = Compliance[2,2]
            S44:float = Compliance[3,3]
            S12:float = Compliance[0,1]
            S13:float = Compliance[0,2]
            BR:float = 1./(2.*(S11+S22+2.*S13)+S33)
            GR:float = 30./((19.*S11-11.*S12-16.*S13)+(8.*S33+12.*S44))
            self.Bmod:float = 0.5*(BV+BR)
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = 0.5*(1.-self.Ymod/(3.*self.Bmod))
            print("* Calculated values by DFT *")
            print("* C11 (GPa): ", self.C11)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C12 (GPa): ", self.C12)
            print("* C13 (GPa): ", self.C13)
            print("* C24 (GPa): ", self.C24)
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
            print("* There are 13 independent elastic constants: ")
            print("* C11, C22, C33, C44, C55, C66, C12, C13, C15, C23, C25, C35, C46")
            """ 0.5 * C11 """
            sym:str = "mono-x"
            dfmat:float = np.array([np.diag([1.+d,1.,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E10:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E1:float = E10[2]
            """ 0.5 * C22 """
            sym:str = "mono-y"
            dfmat:float = np.array([np.diag([1.,1.+d,1.]).tolist() for d in ep])
            emod.calc_Emod(self, ep, dfmat, sym)
            E20:float = emod.Emod_fit(self, para_3rd, sym, sw_3rd=True)
            E2:float = E20[2]
            """ 0.5 * C33 """
            sym:str = "mono-z"
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
            #### you need more distorsions

            self.C11:float = prms.AU2GPa * 2. * E1
            self.C22:float = prms.AU2GPa * 2. * E2
            self.C33:float = prms.AU2GPa * 2. * E3
            self.C44:float = prms.AU2GPa * 2. * E4
            self.C55:float = prms.AU2GPa * 2. * E5
            self.C66:float = prms.AU2GPa * 2. * E6
            self.C12:float = prms.AU2GPa * (E7 - E1 - E2)
            self.C13:float = prms.AU2GPa * (E8 - E1 - E3)
            #self.C15:float =
            #self.C23:float =
            #self.C25:float =
            #self.C35:float =
            #self.C46:float = 
            Econst:float = np.array([self.C11,self.C12,self.C13,self.C22,self.C33,self.C44,self.C55,self.C66])
            BV:float = (self.C11+self.C22+self.C33+2.*(self.C12+self.C13+self.C23))/9.
            GV:float = (self.C11+self.C22+self.C33+3.*(self.C44+self.C55+self.C66)-(self.C12+self.C13+self.C23))/15.
            a:float = self.C33*self.C55 - self.C35**2.
            b:float = self.C23*self.C55 - self.C25*self.C35
            c:float = self.C13*self.C35 - self.C15*self.C33
            d:float = self.C13*self.C55 - self.C15*self.C35
            e:float = self.C13*self.C25 - self.C15*self.C23
            f:float = self.C11*(self.C22*self.C55 - self.C25**2.) - self.C12*(self.C12*self.C55 - self.C15*self.C25) + self.C15*(self.C12*self.C25 - self.C15*self.C22) + self.C25*(self.C23*self.C35 - self.C25*self.C33)
            g:float = self.C11*self.C22*self.C33 - self.C11*self.C23**2. - self.C22*self.C13**2. - self.C33*self.C12**2. + 2.*self.C12*self.C13*self.C23
            Omg:float = 2.*(self.C15*self.C25*(self.C33*self.C12-self.C13*self.C23)+self.C15*self.C35*(self.C22*self.C13-self.C12*self.C23)+self.C25*self.C35
                            *(self.C11*self.C23-self.C12*self.C13))-((self.C15**2.)*(self.C22*self.C33-self.C23**2.)+(self.C25**2.)*(self.C11*self.C33-self.C13**2.)+(self.C35**2.)
                                                                     *(self.C11*self.C22-self.C12**2.)) + g*self.C55
            BR:float = Omg/(a*(self.C11+self.C22-2.*self.C12)+b*(2.*self.C12-2.*self.C11-self.C23)+c*(self.C15-2.*self.C25)+d*(2.*self.C12+2.*self.C23-self.C13-2.*self.C22)+2.*e*(self.C25-self.C15)+f)
            coff:float = 4.*(a*(self.C11+self.C22+self.C12)+b*(self.C11-self.C12-self.C23)+c*(self.C15+self.C25)+d*(self.C22-self.C12-self.C23-self.C13)+e*(self.C15-self.C25)+f)
            GR:float = 15./((coff/Omg) + 3.*((g/Omg) + (self.C44+self.C66)/(self.C44*self.C66-self.C46**2.)))
            self.Bmod:float = 0.5*(BV+BR)
            self.Gmod:float = 0.5*(GV+GR)
            self.Ymod:float = 9.*self.Bmod*self.Gmod/(3.*self.Bmod+self.Gmod)
            self.nu:float = 0.5*(1.-self.Ymod/(3.*self.Bmod))
            print("* Calculated values by DFT *")
            print("* C11 (GPa): ", self.C11)
            print("* C22 (GPa): ", self.C22)
            print("* C33 (GPa): ", self.C33)
            print("* C44 (GPa): ", self.C44)
            print("* C55 (GPa): ", self.C55)
            print("* C66 (GPa): ", self.C66)
            print("* C12 (GPa): ", self.C12)
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

        else:
            print("*** ERROR in emod.get_Econst: brav should be 'cub', 'tet', 'ortho', 'hex', 'trig1', or 'mono'!!!")
            print("*** Other crystal systems (triclinic) are not implemented!!!")
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
    def get_eig(self, path:str = "."):
        """ get eigenvalues from bnd***.spin* in ecalj and EIGENVAL in VASP """
    
        self.eigu:float = [ [] for i in range(prms.nkpath*prms.nkpt) ]
        self.eigd:float = [ [] for i in range(prms.nkpath*prms.nkpt) ]
        self.ehomo:float = -1000.0
        self.elumo:float = 1000.0
    
        if ( prms.code == "ecalj" ):
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
                        print("*** ERROR in cls_emod.get_eig: "+fn+" doesn't exist!!!")
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

        elif ( prms.code == "qe" ):
            print("*** ERROR in cls_emod.get_eig: Sorry, I do not implement 'qe' version now.")
            sys.exit()
            
        else:
            print("*** ERROR in cls_emod.get_eig: code should be 'ecalj', or 'qe'!!!")
            sys.exit()

        # E_g = LUMO - HOMO
        self.egap: float = self.elumo - self.ehomo
        print("* Energy gap (eV): ", self.egap)
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def job_auto(self):
        """ perform jobs of elastic moduli & energy band gaps using ecalj """

        p = pathlib.Path("LISTS")
        if ( not p.exists() ):
            print("*** LISTS does not exist. Skip job_auto.")
            mat_lists:str = []
            brav_lists:str = []
        else:
            mat_lists, brav_lists = np.loadtxt("LISTS",dtype=str,unpack=True,ndmin=0)
            
        for i in range(len(mat_lists)):
            prms.set_prms(mat_lists[i],brav_lists[i])
            sub.run(["mkdir -p RESULTS/"+prms.mat],shell=True)
            sub.run(["cp POS/POSCAR_{mat} RESULTS/{mat}/POSCAR0".format(mat=prms.mat)],shell=True)
            os.chdir("RESULTS/"+prms.mat)
            sub.run(["vasp2ctrl POSCAR0"],shell=True)
            sub.run(["mv ctrls.POSCAR0.vasp2ctrl ctrls."+prms.mat],shell=True)
            (alat,plat,elems,nele,natm,pos,vol) = prms.get_POSCAR("POSCAR0")
            if ( "F" in elems ):
                sub.run(["sed -i -e 's/F/F0/g' ctrls."+prms.mat],shell=True)
                sub.run(["sed -i -e 's/F0e/Fe/g' ctrls."+prms.mat],shell=True)
                with open("ctrls."+prms.mat,"a") as f:
                    f.write("\n")
                    f.write("SPEC \n")
                    for ele in elems:
                        if ( ele == "F" ):
                            f.write("   ATOM=F0 Z=9 \n")
                        else:
                            f.write("   ATOM={ele} Z={Znum} \n".format(ele=ele,Znum=prms.Znum[ele]))
            sub.run(["ctrlgenM1.py "+prms.mat],shell=True)
            sub.run(['mv ctrlgenM1.ctrl.{mat} ctrl.{mat}'.format(mat=prms.mat)], shell=True)
            sub.run(['getsyml '+prms.mat], shell=True)

            sub.run(["{exe}/lmfa {mat} >& llmfa".format(exe=prms.exe,mat=prms.mat)],shell=True)
            sub.run(["mpirun -np {nc} {exe}/lmf {mat} >& llmf".format(nc=prms.nc,exe=prms.exe,mat=prms.mat)],shell=True)
            sub.run(["mkdir LDA"],shell=True)
            sub.run(["{exe}/job_band {mat} -np {nc}".format(exe=prms.exe,mat=prms.mat,nc=prms.nc)],shell=True)
            sub.run(["{exe}/job_pdos {mat} -np {nc}".format(exe=prms.exe,mat=prms.mat,nc=prms.nc)],shell=True)
            sub.run(["cp *.{mat} LDA".format(mat=prms.mat)],shell=True)
            sub.run(["mkGWinput {mat}".format(mat=prms.mat)],shell=True)
            sub.run(["cp GWinput.tmp GWinput"],shell=True)
            sub.run(["{exe}/gwsc {itr} {mat} -np {nc}".format(exe=prms.exe,itr=prms.GWitr,mat=prms.mat,nc=prms.nc)],shell=True)
            
    ### ----------------------------------------------------------------------------- ###
    def get_kmesh(self):
        """ get appropriate kmesh """

        (alat,plat,elems,nele,natm,pos,vol) = prms.get_POSCAR("POSCAR0")
        bvec:float = alat*alat * (2.*np.pi/vol) * np.array([np.cross(plat[1],plat[2]),np.cross(plat[2],plat[0]),np.cross(plat[0],plat[1])])
        bvol:float = abs(np.dot(bvec[0],np.cross(bvec[1],bvec[2])))
        ### standard ==>> Si 12*12*12
        
    ### ----------------------------------------------------------------------------- ###
    def save_emod(self):
        """ save emod data to a file """

        fn:str = "DATA_emod"
        p = pathlib.Path(fn)
        if ( not p.exists() ):
            string:str = "### 1:mat   2:Bmod(GPa)   3:Gmod(GPa)   4:Ymod(GPa)   5:Debye_temp(K)   6:egap(eV)  \n"
            with open(fn,"w") as f:
                f.write(string)
        string:str = "{mat}   {Bmod:.f5}   {Gmod:.f5}   {Ymod:.f5}   {Debye:.f5}   {egap:.f5} \n".format(mat=prms.mat,Bmod=self.Bmod,Gmod=self.Gmod,Ymod=self.Ymod,Debye=self.ThetaD,egap=self.egap)
        with open(fn,"a") as f:
            f.write(string)
