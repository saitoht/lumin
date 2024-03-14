### See "M. A. Reshchikov & H. Morkoc, Luminescence properties of defects in GaN, J. Appl. Phys. 97 (2005) 061301."
###     "Y. Jia, A. Miglio, S. Ponce, M. Mikami, & X. Gonze, First-principles study of the luminescence of Eu2+-doped phosphors, PRB 96 (2017) 125132."
###
### Things to do
### - include anharmonic effects in Multi-D CCD model (maybe by Self-consistent Phonon Theory?)
### - make the code more understandable

import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import os, sys, pathlib
import cls_subs as subs

class ccd:
    """ Class: configuration coordinate diagram """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        const = subs.subs()
        prms = subs.get_prms()

    ### ----------------------------------------------------------------------------- ###
    def run_ccd(self):
        """ execute ccd calculation program """
        
        # For plotting
        plt.rcParams['font.family'] = 'Helvetica'
        plt.rcParams["xtick.labelsize"]=15.0
        plt.rcParams["ytick.labelsize"]=15.0
        plt.rcParams["xtick.major.pad"] = 5
        plt.rcParams["ytick.major.pad"] = 5
        plt.rcParams["axes.labelsize"] = 15.0
        plt.rcParams["axes.linewidth"] = 1.0
        plt.rcParams["axes.labelpad"] = 6
        plt.rcParams["xtick.direction"] = "in" 
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams["xtick.major.width"] = 1.0
        plt.rcParams["ytick.major.width"] = 1.0
        plt.rcParams["xtick.minor.width"] = 0.5
        plt.rcParams["ytick.minor.width"] = 0.5
        plt.rcParams["xtick.major.size"] = 4.5
        plt.rcParams["ytick.major.size"] = 4.5
        plt.rcParams["xtick.minor.size"] = 3.0
        plt.rcParams["ytick.minor.size"] = 3.0
    
        print("* --- START CALCULATIONS OF STOKES SHIFT--- *")
                
        if ( prms.sw_unit in {"eV","nm"} ):
            unit:float = 1.0
        elif ( prms.sw_unit == "cm^-1" ):
            unit:float = const.eV2cm
        get_Stokes()
        get_DeltaQ()
        if ( sw_eg ):
            calc_Eeg(prms.stateg)
            calc_Eeg(prms.statee)
            (Qming_3dim,Omegag,Sem,Qming_2dim,Omegag_2dim,Sem_2dim) = fit_Ecurve(prms.stateg,prms.EFCg)
            (Qmine_3dim,Omegae,Sabs,Qmine_dim,Omegae_2dim,Sabs_2dim) = fit_Ecurve(prms.statee,self.EFCe)
        
        Line_shape()
        Ecenter_shift()
        FWHM()
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            print("* Absorption energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eabs0*unit)
            print("* Emission energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eem0*unit)
            print("* EFCg = Ee - Eg ({sunit}): ".format(sunit=prms.sw_unit), prms.EFCg*unit)
            print("* EFCe = Eg* - Ee* ({sunit}): ".format(sunit=prms.sw_unit), self.EFCe*unit)
            print("* Stokes shift (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.deltaS*unit)
            print("* Zero phonon photoemission line ({sunit}): ".format(sunit=prms.sw_unit), self.EZPL*unit)
            print("* FWHM (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.W0*unit)
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3]: ", self.dQvec)
            print("* hbar*Omegag ({sunit}): ".format(sunit=prms.sw_unit), self.Omegag*unit)
            print("* hbar*Omegae ({sunit}): ".format(sunit=prms.sw_unit), self.Omegae*unit)
        elif ( prms.sw_unit == "nm" ):
            print("* Absorption energy (nm): ", subs.E2lambda(prms.Eabs0))
            print("* Emission energy (nm): ", subs.E2lambda(prms.Eem0))
            print("* EFCg = Ee - Eg (nm): ", subs.E2lambda(prms.EFCg))
            print("* EFCe = Eg* - Ee* (nm): ", subs.E2lambda(self.EFCe))
            print("* Stokes shift (@ 0K) (nm): ", subs.E2lambda(self.deltaS))
            print("* Zero phonon photoemission line (nm): ", subs.E2lambda(self.EZPL))
            print("* FWHM (@ 0K) (nm): ", subs.E2lambda(self.W0))
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3]: ", self.dQvec)
            print("* hbar*Omegag (nm): ", subs.E2lambda(self.Omegag))
            print("* hbar*Omegae (nm): ", subs.E2lambda(self.Omegae))
        else:
            print("* ERROR: sw_unit should be 'eV', 'nm', or 'cm^-1'!!!")
            sys.exit()
    
        print("* Sabs: ", self.Sabs)
        print("* Sem: ", self.Sem)
        print("*")
    
        if ( prms.sw_plt_ccd ):
            print("* --- PLOT LINE SHAPE --- *")
            size:int = 30
            plt.ylabel("Intensity (a.u.)")
            plt.axhline(0, lw=0.5, c="black", linestyle="dashed")
            if ( prms.sw_unit in {"eV","cm^-1"} ):
                plt.xlabel("Energy ({sunit})".format(sunit=prms.sw_unit))
                plt.plot(unit*energy, self.Labs, lw=1.0, c="black")
                plt.scatter(prms.Eabs0*unit,0.0,c="black",s=size,marker="o")
                plt.plot(unit*energy, self.Lem, lw=1.0, c="red")
                plt.scatter(prms.Eem0*unit,0.0,c="red",s=size,marker="o")
            elif ( prms.sw_unit == "nm" ):
                plt.xlabel("Wave length (nm)")
                plt.plot(subs.E2lambda(energy), self.Labs, lw=1.0, c="black")
                plt.plot(subs.E2lambda(energy), self.Lem, lw=1.0, c="red")
            plt.savefig("Spectrum.pdf")
            plt.show()
        
            print("* --- PLOT TEMPERATURE DEPENDENCE OF E_EM & E_ABS --- *")
            plt.xlabel("Temperature (K)")
            plt.ylabel("Emission energy ({sunit})".format(sunit=prms.sw_unit))
            if ( prms.sw_unit in {"eV","cm^-1"} ):
                plt.plot(self.temp, unit*self.Eabs, lw=1.0, c="black")
                plt.plot(self.temp, unit*self.Eem, lw=1.0, c="red")
            elif ( prms.sw_unit == "nm" ):
                plt.plot(self.temp, subs.E2lambda(self.Eabs), lw=1.0, c="black")
                plt.plot(self.temp, subs.E2lambda(self.Eem), lw=1.0, c="red")
            plt.savefig("Epeak_Temp.pdf")
            plt.show()
        
            print("* --- PLOT TEMPERATURE DEPENDENCE OF STOKES SHIFT --- *")
            plt.xlabel("Temperature (K)")
            plt.ylabel("Stokes shift ({sunit})".format(sunit=prms.sw_unit))
            if ( prms.sw_unit in {"eV","cm^-1"} ):
                plt.plot(self.temp, unit*(self.Eabs-self.Eem), lw=1.0, c="black")
            elif ( prms.sw_unit == "nm" ):
                plt.plot(self.temp, subs.E2lambda(self.Eabs)-subs.E2lambda(self.Eem), lw=1.0, c="black")
            plt.savefig("Stokes_Temp.pdf")
            plt.show()
        
            print("* --- PLOT TEMPERATURE DEPENDENCE OF FWHM --- *")
            plt.xlabel("Temperature (K)")
            plt.ylabel("FWHM ({sunit})".format(sunit=prms.sw_unit))
            plt.plot(self.temp, unit*self.W, lw=1.0, c="black")
            plt.savefig("FWHM_Temp.pdf")
            plt.show()
        else:
            print("* --- NO PLOT --- *")
        
        print("* --- FINISH CALCULATIONS OF STOKES SHIFT--- *")

    ### ----------------------------------------------------------------------------- ###
    def calc_Eeg(self, state:str):
        """ calculate total enegies of intermediate states """
        
        print("* --- gernerate POSCAR in between ground & excited states --- *")
        (alat_g, plat_g, elements_g, nelems_g, natm_g, pos_g) = get_POSCAR("POSCAR_"+prms.stateg)
        (alat_e, plat_e, elements_e, nelems_e, natm_e, pos_e) = get_POSCAR("POSCAR_"+prms.statee)
        Qarr: float = np.linspace(-prms.dQ, 1.+prms.dQ, prms.ndiv_eg)
        for i, Q in enumerate(Qarr):
            Rmid: float = np.zeros((len(pos_g),3))
            Rmid[:,:] = (Q / self.deltaQ) * (pos_e[:,:] - pos_g[:,:]) + pos_g[:,:]
            fn_eg: str = "POSCAR_eg{idn}".format(idn=i+1)
            sub.run(["head -8 POSCAR_g > {fn}".format(fn=fn_eg)], shell=True)
            strcoord: str = ""
            for j in range(len(Rmid)):
                strcoord += " {pos_x}  {pos_y}  {pos_z} \n".format(pos_x=Rmid[j,0], pos_y=Rmid[j,1], pos_z=Rmid[j,2])
            with open(fn_eg, "a") as f:
                f.write(strcoord)
            print("* Q={Q} {fn} created!".format(Q=Q, fn=fn_eg))
        print("* Finish!")
        print("*")
    
        fn: str = "DATA_"+state+".dat"
        fp = pathlib.Path(fn)
        if ( not fp.exists() ):
            string: str = "### Q((amu)^(1/2)*ang)   Etot(eV) \n"
            with open(fn,"w") as f:
                f.write(string)
        for i, qc in enumerate(Qarr):
            ndir: str = "eg{num}".format(num=i+1)
            pd = pathlib.Path(ndir+"/"+state)
            if ( not pd.is_dir() ):
                sub.run(["mkdir -p "+ndir+"/"+state], shell=True)
                os.chdir(ndir+"/"+state)
                sub.run(["cp ../../header_"+state+".in "+prms.mat+"-eg"+str(i+1)+".scf.in"], shell=True)
                (alat, plat, elements, nelems, natm, pos) = subs.get_POSCAR("../../POSCAR_eg"+str(i+1))
                ic: int = 0
                string: str = ""
                for j, nele in enumerate(nelems):
                    for k in range(nele):
                        string += " "+elements[j]+"  {x}  {y}  {z}  \n".format(x=pos[ic,0],y=pos[ic,1],z=pos[ic,2])
                        ic += 1
                with open(mat+"-eg"+str(i+1)+".scf.in","a") as f:
                    f.write(string)
                # run scf calculation
                sub.run(["mpirun -np "+str(prms.nc)+" "+prms.exe+"/pw.x < "+prms.mat+"-eg"+str(i+1)+".scf.in > "+prms.mat+"-eg"+str(i+1)+".scf.out"], shell=True)
                sub.run(["rm work/"+prms.mat+".save/wfc*"], shell=True)
                sub.run(["grep ! "+prms.mat+"-eg"+str(i+1)+".scf.out > grep.out"], shell=True)
                data: str = np.loadtxt("grep.out",dtype="str",unpack=True,ndmin=0)
                os.chdir("../../")
                string = "{Q}   {Etot} \n".format(Q=qc, Etot=data[4])
                with open(fn, "a") as f:
                    f.write(string)

    ### ----------------------------------------------------------------------------- ###
    def fit_Ecurve(self, state:str, EFC:float):
        """ fiting Energy curve by polynomial """
        
        print("* --- Fitting energy curve --- *")
        fn:str = "DATA_"+stat+".dat"
        Q, Etot = np.loadtxt(fn,dtype='float',unpack=True,ndmin=0)
        x:float = np.linspace(-0.1+Q[0],0.1+Q[len(Q)-1],1000)
        weight:float = np.ones(len(Q))
        coef_3dim:float = np.polyfit(Q, (Etot-min(Etot))*const.Ry, 3, w=weight)
        Efit_3dim:float = np.poly1d(coef_3dim)(x)
        coef_2dim:float = np.polyfit(Q, (Etot-min(Etot))*const.Ry, 2, w=weight)
        Efit_2dim:float = np.poly1d(coef_2dim)(x)

        # consider Anharmonic effects
        cfunit:float = const.hbar*1.e10*np.sqrt(1./(const.ep*const.uatm))
        Qmin_3dim:float = (-coef_3dim[1]+np.sqrt(coef_3dim[1]**2.-3.*coef_3dim[0]*coef_3dim[2]))/(3.*coef_3dim[0])
        Omega_3dim:float = cfunit*np.sqrt(6.*coef_3dim[0]*Qmin_3dim+2.*coef_3dim[1])
        S_3dim:float = EFC / Omega_3dim
        print("* Fitted by cubic equation: ")
        print(coef_3dim)
        print("* hbar*Omega (eV) [cubic]: sqrt({cube}Q+{square})".format(cube=(cfunit**2.)*6.*coef_3dim[0], square=(cfunit**2.)*2.*coef_3dim[1]))
        print("* Qmin (amu^1/2 ang): {Qmin}".format(Qmin=Qmin_3dim))
        print("* hbar*Omega (eV) @ Q=Qmin: {cubic}".format(cubic=Omega_3dim))
        print("*")
    
        # Harmonic approximation
        Omega_2dim:float = cfunit*np.sqrt(2.*coef_2dim[0])
        Qmin_2dim:float = - coef_2dim[1] / (2.*coef_2dim[0])
        S_2dim:float = EFC / Omega_2dim
        print("* Fitted by quadratic equation: ")
        print(coef_2dim)
        print("* hbar*Omega (eV) [quadratic]: {square}".format(square=Omega_2dim))
        print("* Qmin (amu^1/2 ang): {Qmin}".format(Qmin=Qmin_2dim))
    
        plt.xlabel(r"Q ($\sqrt{\mathrm{amu}} \cdot \AA$)")
        plt.ylabel("Energy (eV)")
        plt.scatter(Q,(Etot-min(Etot))*const.Ry,marker="o",label="data",edgecolor="darkred",color="white",s=30)
        plt.plot(x,Efit_2dim,linestyle="dashed",color="mediumblue",label="fit 2dim")
        plt.scatter(Qmin_2dim, np.poly1d(coef_2dim)(Qmin_2dim), marker="*", color="darkblue", s=60)
        plt.plot(x,Efit_3dim,linestyle="dashed",color="coral",label="fit 3dim")
        plt.scatter(Qmin_3dim, np.poly1d(coef_3dim)(Qmin_3dim), marker="*", color="crimson", s=60)
        plt.legend()
        plt.savefig("Ecurve_fit_"+state+".pdf")
        if ( prms.sw_plt_ccd ):
            plt.show()
        return (Qmin_3dim, Omega_3dim, S_3dim, Qmin_2dim, Omega_2dim, S_2dim)
        
    ### ----------------------------------------------------------------------------- ###
    def get_Stokes(self):
        """ Stokes shift, Franck-Condon parameter, and zero-phonon energy """
        
        self.deltaS: float = prms.Eabs0 - prms.Eem0
        self.EFCe: float = self.deltaS - prms.EFCg
        self.EZPL: float = prms.Eem0 + prms.EFCg

    ### ----------------------------------------------------------------------------- ###
    def get_DeltaQ(self):
        """ normal coordinate difference DetaQ """
        
        pg = pathlib.Path("POSCAR_"+prms.stateg)
        pe = pathlib.Path("POSCAR_"+prms.statee)
        if ( not pg.exists() ):
            print("*** ERROR: POSCAR_g doesn't exist!!!")
            sys.exit()
        if ( not pe.exists() ):
            print("*** ERROR: POSCAR_e doesn't exist!!!")
            sys.exit()
        (alat_g, plat_g, elements_g, nelems_g, natm_g, pos_g, volume) = subs.get_POSCAR("POSCAR_"+prms.stateg)
        (alat_e, plat_e, elements_e, nelems_e, natm_e, pos_e, volume) = subs.get_POSCAR("POSCAR_"+prms.statee)
        if ( not abs(alat_g - alat_e) < 1.0e-5 ):
            print("*** ERROR: different alat between the grond state and the excited state!!!")
            sys.exit()
            
        mass:float = [const.ELEMS_MASS[ele]*const.uatm/const.me for ele in elements_g]
        alat_g = alat_g / const.Bohr
        alat_e = alat_e / const.Bohr
        pos_g = alat_g * np.dot( pos_g, plat_g )
        pos_e = alat_e * np.dot( pos_e, plat_e )
        deltaQ_sq:float = 0.0
        deltaR_sq:float = 0.0
        dQvec_sq:float = np.zeros(3)
        count:int = 0
        for i, nel in enumerate(nelems_g):
            for j in range(nel):
                deltaQ_sq += mass[i] * ( sum([(pos_e[count,k]-pos_g[count,k])**2. for k in range(3)]) )
                deltaR_sq += sum([(pos_e[count,k]-pos_g[count,k])**2. for k in range(3)])
                dQvec_sq[:] += mass[i] * (pos_e[count,:]-pos_g[count,:])**2.
                count += 1
        self.deltaQ:float = np.sqrt(const.me/const.uatm)*const.Bohr * np.sqrt(deltaQ_sq)
        self.dQvec:float = np.sqrt(const.me/const.uatm)*const.Bohr * np.array([np.sqrt(dQvec_sq[i]) for i in range(3)])
        self.deltaR:float = const.Bohr * np.sqrt(deltaR_sq)
        self.M:float = deltaQ_sq / deltaR_sq

        cfunit:float = const.hbar*1.e10*np.sqrt(1./(const.ep*const.uatm))
        self.Omegag:float = cfunit * np.sqrt(2. * prms.EFCg / (self.deltaQ**2.)) 
        self.Omegae:float = cfunit * np.sqrt(2. * self.EFCe / (self.deltaQ**2.)) 
        self.Sabs:float = self.EFCe / self.Omegae
        self.Sem:float = prms.EFCg / self.Omegag
        
    ### ----------------------------------------------------------------------------- ###
    def Line_shape(self):
        """ Line shape of optical spectra """
        
        if ( prms.sw_unit in {"eV","nm"} ):
            self.energy:float = np.linspace(prms.emin_plt, prms.emax_plt, prms.ndive)
        else:
            self.energy:float = np.linspace(subs.E2lambda(prms.emax_plt),subs.E2lambda(prms.emin_plt),prms.ndive)
            self.energy = subs.lambda2E(self.energy)
        self.Lem:float = np.zeros(prms.ndive)
        self.Labs:float = np.zeros(prms.ndive)
        for n in range(prms.nmax):
            self.Lem[:] += ( np.exp(-self.Sem)*self.Sem**float(n) / float(np.math.factorial(n)) ) * subs.Lorentzian(self.EZPL - float(n)*self.Omegag - self.energy[:])
            self.Labs[:] += ( np.exp(-self.Sabs)*self.Sabs**float(n) / float(np.math.factorial(n)) ) * subs.Lorentzian(self.EZPL + float(n)*self.Omegae - self.energy[:])
        self.Lem = self.I0 * self.Lem / max(self.Lem)
        self.Labs = self.I0 * self.Labs / max(self.Labs)

    ### ----------------------------------------------------------------------------- ###
    def Ecenter_shift(self)
        """ Eem shift and Eabs shift as a function of temperature """
        
        ### revise to prms.tempmin in the future
        self.tempmin:float = 1.
        self.tempmax:float = 1.e3
        self.temp:float = np.linspace(self.tempmin, self.tempmax, prms.ndivtemp)
        self.Eabs:float = prms.Eabs0 + ((self.Omegae**2. - self.Omegag**2.)/(self.Omegag**2.)) * const.kB*self.temp
        self.Eem:float = prms.Eem0 + ((self.Omegag**2. - self.Omegae**2.)/(self.Omegae**2.) +
                                       (8.*(self.Omegag**4.)*self.deltaS)/(self.Omegae**2.*(self.Omegag**2.+self.Omegae**2.)*prms.Eem0)) * prms.kB*self.temp

    ### ----------------------------------------------------------------------------- ###
    def FWHM(self):
        """ Full width half maximum of transitions and its temperature dependence """
        
        self.temp:float = np.linspace(self.tempmin, self.tempmax, prms.ndivtemp)
        self.W0:float = self.Sem * self.Omegag * np.sqrt(8.0*np.log(2.)) / np.sqrt(self.Sabs)
        if ( prms.sw_unit == "nm" ):
            W0min:float = prms.Eem0 - 0.5*self.W0
            W0max:float = prms.Eem0 + 0.5*self.W0
            self.W0 = subs.E2lambda(W0min) - subs.E2lambda(W0max)
        self.W:float = self.W0 * np.sqrt( 1. / np.tanh(self.Omegae/(2.*const.kB*self.temp)) )
