import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import os, sys, pathlib
import cls_subs as subs

prms = subs.subs()

class ccd:
    """ Class: configuration coordinate diagram """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor of configuration coordinate model """
        
        print("* --- Start calculation of configuration coordinate model--- *")
        if ( prms.sw_unit in {"eV","nm"} ):
            self.unit:float = 1.0
        elif ( prms.sw_unit == "cm^-1" ):
            self.unit:float = prms.eV2cm
        ccd.get_Stokes(self)
        ccd.get_DeltaQ(self)
        if ( prms.sw_eg ):
            ccd.calc_Eeg(self,prms.stateg)
            ccd.calc_Eeg(self,prms.statee)
            (Qming_3dim,Omegag,Sem,Qming_2dim,Omegag_2dim,Sem_2dim) = ccd.fit_Ecurve(self,prms.stateg,prms.EFCg)
            (Qmine_3dim,Omegae,Sabs,Qmine_dim,Omegae_2dim,Sabs_2dim) = ccd.fit_Ecurve(self,prms.statee,self.EFCe)
        ccd.calc_Line_shape(self)
        ccd.calc_Ecenter_shift(self)
        ccd.calc_FWHM(self)

        print("*** check the parameters in ccd ***")
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            print("* Absorption energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eabs0*self.unit)
            print("* Emission energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eem0*self.unit)
            print("* EFCg = Ee - Eg ({sunit}): ".format(sunit=prms.sw_unit), prms.EFCg*self.unit)
            print("* EFCe = Eg* - Ee* ({sunit}): ".format(sunit=prms.sw_unit), self.EFCe*self.unit)
            print("* Stokes shift (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.deltaS*self.unit)
            print("* Zero phonon photoemission line ({sunit}): ".format(sunit=prms.sw_unit), self.EZPL*self.unit)
            print("* FWHM (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.W0*self.unit)
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3]: ", self.dQvec)
            print("* hbar*Omegag ({sunit}): ".format(sunit=prms.sw_unit), self.Omegag*self.unit)
            print("* hbar*Omegae ({sunit}): ".format(sunit=prms.sw_unit), self.Omegae*self.unit)
        elif ( prms.sw_unit == "nm" ):
            print("* Absorption energy (nm): ", prms.E2lambda(prms.Eabs0))
            print("* Emission energy (nm): ", prms.E2lambda(prms.Eem0))
            print("* EFCg = Ee - Eg (nm): ", prms.E2lambda(prms.EFCg))
            print("* EFCe = Eg* - Ee* (nm): ", prms.E2lambda(self.EFCe))
            print("* Stokes shift (@ 0K) (nm): ", prms.E2lambda(self.deltaS))
            print("* Zero phonon photoemission line (nm): ", prms.E2lambda(self.EZPL))
            print("* FWHM (@ 0K) (nm): ", prms.E2lambda(self.W0))
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3]: ", self.dQvec)
            print("* hbar*Omegag (nm): ", prms.E2lambda(self.Omegag))
            print("* hbar*Omegae (nm): ", prms.E2lambda(self.Omegae))
        else:
            print("* ERROR: sw_unit should be 'eV', 'nm', or 'cm^-1'!!!")
            sys.exit()
        print("* Sabs: ", self.Sabs)
        print("* Sem: ", self.Sem)
        print("*")
    
        if ( prms.sw_plt_ccd ):
            ccd.plt_ccd(self)
        else:
            print("* --- NO PLOT --- *")
        print("* --- FINISH CALCULATIONS OF STOKES SHIFT--- *")

    ### ----------------------------------------------------------------------------- ###
    def plt_ccd(self):
        """ plot the results of configuration coordinate model """

        def Eg(x,x0,y0):
            return (y0/(x0**2.)) * (x**2.)

        def Ee(x,x0,y0,y1):
            return ((y1-y0)/(x0**2.)) * ((x-x0)**2.) + y0

        plt.rcParams['font.family'] = 'Helvetica'
        plt.rcParams["xtick.labelsize"] = 15.0
        plt.rcParams["ytick.labelsize"] = 15.0
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

        print("*** PLOT 1D-CCD ***")
        xQ:float = np.linspace(-0.5*self.deltaQ,1.5*self.deltaQ,1000)
        yEg:float = Eg(xQ,self.deltaQ,prms.EFCg)
        yEe:float = Ee(xQ,self.deltaQ,prms.EFCg+prms.Eem0,prms.Eabs0)
        Qarr:float = [0.0, 0.0, self.deltaQ, self.deltaQ]
        Earr:float = [0.0, prms.Eabs0, prms.EFCg, prms.EFCg+prms.Eem0]
        plt.figure(figsize=(6,6.5))
        plt.xlim(-0.5*self.deltaQ,1.5*self.deltaQ)
        plt.ylim(-0.1, prms.emax_ccd)
        plt.xlabel(r"$\Delta Q$ (amu$^{1/2}\cdot\mathrm{\AA}$)")
        plt.ylabel("Energy (eV)")
        plt.scatter(Qarr, Earr, color="white", marker="o", edgecolor="mediumblue", s=80)
        plt.plot(xQ,yEg,color="red",lw=1.0)
        plt.plot(xQ,yEe,color="red",lw=1.0)
        plt.plot([0.0,0.0],[-0.1,prms.Eabs0],color="black",linestyle="dotted",lw=0.5)
        plt.plot([self.deltaQ,self.deltaQ],[-0.1,prms.EFCg+prms.Eem0],color="black",linestyle="dotted",lw=0.5)
        plt.plot([-0.5*self.deltaQ,self.deltaQ],[prms.EFCg+prms.Eem0,prms.EFCg+prms.Eem0],color="black",linestyle="dotted",lw=0.5)
        plt.plot([-0.5*self.deltaQ,self.deltaQ],[prms.EFCg,prms.EFCg],color="black",linestyle="dotted",lw=0.5)
        plt.plot([-0.5*self.deltaQ,1.5*self.deltaQ],[0.0,0.0],color="black",linestyle="dotted",lw=0.5)
        plt.plot([-0.5*self.deltaQ,0.0],[prms.Eabs0,prms.Eabs0],color="black",linestyle="dotted",lw=0.5)
        plt.savefig("1DCCD_{state}.pdf".format(state=prms.statee))
        plt.show()
        
        print("*** PLOT LINE SHAPE ***")
        size:int = 30
        plt.ylabel("Intensity (a.u.)")
        plt.axhline(0, lw=0.5, c="black", linestyle="dashed")
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            plt.xlabel("Energy ({sunit})".format(sunit=prms.sw_unit))
            plt.plot(self.unit*self.energy, self.Labs, lw=1.0, c="black")
            plt.scatter(prms.Eabs0*self.unit,0.0,c="black",s=size,marker="o")
            plt.plot(self.unit*self.energy, self.Lem, lw=1.0, c="red")
            plt.scatter(prms.Eem0*self.unit,0.0,c="red",s=size,marker="o")
        elif ( prms.sw_unit == "nm" ):
            plt.xlabel("Wave length (nm)")
            plt.plot(prms.E2lambda(self.energy), self.Labs, lw=1.0, c="black")
            plt.plot(prms.E2lambda(self.energy), self.Lem, lw=1.0, c="red")
        plt.savefig("Spectrum.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF E_EM & E_ABS ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Emission energy ({sunit})".format(sunit=prms.sw_unit))
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            plt.plot(self.temp, self.unit*self.Eabs, lw=1.0, c="black")
            plt.plot(self.temp, self.unit*self.Eem, lw=1.0, c="red")
        elif ( prms.sw_unit == "nm" ):
            plt.plot(self.temp, prms.E2lambda(self.Eabs), lw=1.0, c="black")
            plt.plot(self.temp, prms.E2lambda(self.Eem), lw=1.0, c="red")
        plt.savefig("Epeak_Temp.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF STOKES SHIFT ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Stokes shift ({sunit})".format(sunit=prms.sw_unit))
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            plt.plot(self.temp, self.unit*(self.Eabs-self.Eem), lw=1.0, c="black")
        elif ( prms.sw_unit == "nm" ):
            plt.plot(self.temp, prms.E2lambda(self.Eabs)-prms.E2lambda(self.Eem), lw=1.0, c="black")
        plt.savefig("Stokes_Temp.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF FWHM ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("FWHM ({sunit})".format(sunit=prms.sw_unit))
        plt.plot(self.temp, self.unit*self.W, lw=1.0, c="black")
        plt.savefig("FWHM_Temp.pdf")
        plt.show()
            
    ### ----------------------------------------------------------------------------- ###
    def calc_Eeg(self, state:str):
        """ calculate total enegies of intermediate states """
        
        print("* --- gernerate POSCAR in between ground & excited states --- *")
        (alat_g,plat_g,elements_g,nelems_g,natm_g,pos_g) = prms.get_POSCAR("POSCAR_"+prms.stateg)
        (alat_e,plat_e,elements_e,nelems_e,natm_e,pos_e) = prms.get_POSCAR("POSCAR_"+prms.statee)
        Qarr:float = np.linspace(-prms.dQ, 1.+prms.dQ, prms.ndiv_eg)
        for i, Q in enumerate(Qarr):
            Rmid:float = np.zeros((len(pos_g),3))
            Rmid[:,:] = (Q / self.deltaQ) * (pos_e[:,:] - pos_g[:,:]) + pos_g[:,:]
            fn_eg: str = "POSCAR_eg{ic}".format(ic=i+1)
            sub.run(["head -8 POSCAR_{state} > {fn}".format(state=prms.stateg, fn=fn_eg)], shell=True)
            strcoord:str = ""
            for j in range(len(Rmid)):
                strcoord += " {pos_x}  {pos_y}  {pos_z} \n".format(pos_x=Rmid[j,0], pos_y=Rmid[j,1], pos_z=Rmid[j,2])
            with open(fn_eg, "a") as f:
                f.write(strcoord)
            print("* Q={Q} {fn} created!".format(Q=Q, fn=fn_eg))
        print("* Finish!")
        print("*")
    
        fn: str = "DATA_{state}.dat".format(state=state)
        fp = pathlib.Path(fn)
        if ( not fp.exists() ):
            string: str = "### Q((amu)^(1/2)*ang)   Etot(eV) \n"
            with open(fn,"w") as f:
                f.write(string)
        for i, qc in enumerate(Qarr):
            ndir: str = "eg{ic}".format(ic=i+1)
            pd = pathlib.Path("{ndir}/{state}".format(ndir=ndir, state=state))
            if ( not pd.is_dir() ):
                sub.run(["mkdir -p {ndir}/{state}".format(ndir=ndir, state=state)], shell=True)
                os.chdir("{ndir}/{state}".format(ndir=ndir, state=state))
                sub.run(["cp ../../header_{state}.in {mat}-eg{ic}.scf.in".format(state=state, mat=prms.mat, ic=i+1)], shell=True)
                (alat, plat, elements, nelems, natm, pos) = prms.get_POSCAR("../../POSCAR_eg{ic}".format(ic=i+1))
                ic: int = 0
                string: str = ""
                for j, nele in enumerate(nelems):
                    for k in range(nele):
                        string += " {ele}  {x}  {y}  {z}  \n".format(ele=elements[j],x=pos[ic,0],y=pos[ic,1],z=pos[ic,2])
                        ic += 1
                with open("{mat}-eg{ic}.scf.in".format(mat=prms.mat,ic=i+1),"a") as f:
                    f.write(string)
                """ run scf calculation """
                sub.run(["mpirun -np {nc} {exe}/pw.x < {mat}-eg{ic}.scf.in > {mat}-eg{ic}.scf.out".format(nc=prms.nc, exe=prms.exe, ic=i+1, mat=prms.mat)], shell=True)
                sub.run(["rm work/{mat}.save/wfc*".format(mat=prms.mat)], shell=True)
                sub.run(["grep ! {mat}-eg{ic}.scf.out > grep.out".format(mat=prms.mat, ic=i+1)], shell=True)
                data: str = np.loadtxt("grep.out",dtype="str",unpack=True,ndmin=0)
                os.chdir("../../")
                string = "{Q}   {Etot} \n".format(Q=qc, Etot=data[4])
                with open(fn, "a") as f:
                    f.write(string)

    ### ----------------------------------------------------------------------------- ###
    def fit_Ecurve(self, state:str, EFC:float):
        """ fiting Energy curve by polynomial """
        
        print("* --- Fitting energy curve --- *")
        fn:str = "DATA_{state}.dat".format(state=state)
        Q, Etot = np.loadtxt(fn,dtype='float',unpack=True,ndmin=0)
        x:float = np.linspace(-0.1+Q[0],0.1+Q[len(Q)-1],1000)
        weight:float = np.ones(len(Q))
        coef_3dim:float = np.polyfit(Q, (Etot-min(Etot))*prms.Ry, 3, w=weight)
        Efit_3dim:float = np.poly1d(coef_3dim)(x)
        coef_2dim:float = np.polyfit(Q, (Etot-min(Etot))*prms.Ry, 2, w=weight)
        Efit_2dim:float = np.poly1d(coef_2dim)(x)

        """ consider Anharmonic effects """
        cfunit:float = prms.hbar*1.e10*np.sqrt(1./(prms.ep*prms.uatm))
        Qmin_3dim:float = (-coef_3dim[1]+np.sqrt(coef_3dim[1]**2.-3.*coef_3dim[0]*coef_3dim[2]))/(3.*coef_3dim[0])
        Omega_3dim:float = cfunit*np.sqrt(6.*coef_3dim[0]*Qmin_3dim+2.*coef_3dim[1])
        S_3dim:float = EFC / Omega_3dim
        print("* Fitted by cubic equation: ")
        print(coef_3dim)
        print("* hbar*Omega (eV) [cubic]: sqrt({cube}Q+{square})".format(cube=(cfunit**2.)*6.*coef_3dim[0], square=(cfunit**2.)*2.*coef_3dim[1]))
        print("* Qmin (amu^1/2 ang): {Qmin}".format(Qmin=Qmin_3dim))
        print("* hbar*Omega (eV) @ Q=Qmin: {cubic}".format(cubic=Omega_3dim))
        print("*")
    
        """ Harmonic approximation """
        Omega_2dim:float = cfunit*np.sqrt(2.*coef_2dim[0])
        Qmin_2dim:float = - coef_2dim[1] / (2.*coef_2dim[0])
        S_2dim:float = EFC / Omega_2dim
        print("* Fitted by quadratic equation: ")
        print(coef_2dim)
        print("* hbar*Omega (eV) [quadratic]: {square}".format(square=Omega_2dim))
        print("* Qmin (amu^1/2 ang): {Qmin}".format(Qmin=Qmin_2dim))
    
        plt.xlabel(r"Q ($\sqrt{\mathrm{amu}} \cdot \AA$)")
        plt.ylabel("Energy (eV)")
        plt.scatter(Q,(Etot-min(Etot))*prms.Ry,marker="o",label="data",edgecolor="darkred",color="white",s=30)
        plt.plot(x,Efit_2dim,linestyle="dashed",color="mediumblue",label="fit 2dim")
        plt.scatter(Qmin_2dim, np.poly1d(coef_2dim)(Qmin_2dim), marker="*", color="darkblue", s=60)
        plt.plot(x,Efit_3dim,linestyle="dashed",color="coral",label="fit 3dim")
        plt.scatter(Qmin_3dim, np.poly1d(coef_3dim)(Qmin_3dim), marker="*", color="crimson", s=60)
        plt.legend()
        plt.savefig("Ecurve_fit_{state}.pdf".format(state=state))
        if ( prms.sw_plt_ccd ):
            plt.show()
        return (Qmin_3dim, Omega_3dim, S_3dim, Qmin_2dim, Omega_2dim, S_2dim)
        
    ### ----------------------------------------------------------------------------- ###
    def get_Stokes(self):
        """ Stokes shift, Franck-Condon parameter, and zero-phonon energy """
        
        self.deltaS:float = prms.Eabs0 - prms.Eem0
        self.EFCe:float = self.deltaS - prms.EFCg
        self.EZPL:float = prms.Eem0 + prms.EFCg

    ### ----------------------------------------------------------------------------- ###
    def get_DeltaQ(self):
        """ normal coordinate difference DetaQ """
        
        pg = pathlib.Path("POSCAR_"+prms.stateg)
        pe = pathlib.Path("POSCAR_"+prms.statee)
        if ( not pg.exists() ):
            print("*** ERROR in ccd.get_DeltaQ: POSCAR_{sym} doesn't exist!!!".format(sym=prms.stateg))
            sys.exit()
        if ( not pe.exists() ):
            print("*** ERROR in ccd.get_DeltaQ: POSCAR_{sym} doesn't exist!!!".format(sym=prms.statee))
            sys.exit()
        (alat_g, plat_g, elements_g, nelems_g, natm_g, pos_g, volume) = prms.get_POSCAR("POSCAR_"+prms.stateg)
        (alat_e, plat_e, elements_e, nelems_e, natm_e, pos_e, volume) = prms.get_POSCAR("POSCAR_"+prms.statee)
        if ( not abs(alat_g - alat_e) < 1.0e-5 ):
            print("*** ERROR in ccd.get_DeltaQ: different alat between the grond state and the excited state!!!")
            sys.exit()
            
        mass:float = [prms.ELEMS_MASS[ele]*prms.uatm/prms.me for ele in elements_g]
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
        self.deltaQ:float = np.sqrt(prms.me/prms.uatm) * prms.Bohr * np.sqrt(deltaQ_sq)
        self.dQvec:float = np.sqrt(prms.me/prms.uatm) * prms.Bohr * np.array([np.sqrt(dQvec_sq[i]) for i in range(3)])
        self.deltaR:float = prms.Bohr * np.sqrt(deltaR_sq)
        self.M:float = deltaQ_sq / deltaR_sq

        cfunit:float = prms.hbar*1.e10*np.sqrt(1./(prms.ep*prms.uatm))
        self.Omegag:float = cfunit * np.sqrt(2. * prms.EFCg / (self.deltaQ**2.)) 
        self.Omegae:float = cfunit * np.sqrt(2. * self.EFCe / (self.deltaQ**2.)) 
        self.Sabs:float = self.EFCe / self.Omegae
        self.Sem:float = prms.EFCg / self.Omegag
        
    ### ----------------------------------------------------------------------------- ###
    def calc_Line_shape(self):
        """ Line shape of optical spectra """
        
        if ( prms.sw_unit in {"eV","nm"} ):
            self.energy:float = np.linspace(prms.emin_ccd, prms.emax_ccd, prms.ndiv_e)
        else:
            self.energy:float = np.linspace(prms.E2lambda(prms.emax_ccd),prms.E2lambda(prms.emin_ccd),prms.ndiv_e)
            self.energy = prms.lambda2E(self.energy)
        self.Lem:float = np.zeros(prms.ndiv_e)
        self.Labs:float = np.zeros(prms.ndiv_e)
        for n in range(prms.nmax):
            self.Lem[:] += ( np.exp(-self.Sem)*self.Sem**float(n) / float(np.math.factorial(n)) ) * prms.Lorentzian(self.EZPL - float(n)*self.Omegag - self.energy[:])
            self.Labs[:] += ( np.exp(-self.Sabs)*self.Sabs**float(n) / float(np.math.factorial(n)) ) * prms.Lorentzian(self.EZPL + float(n)*self.Omegae - self.energy[:])
        self.Lem = prms.I0 * self.Lem / max(self.Lem)
        self.Labs = prms.I0 * self.Labs / max(self.Labs)

    ### ----------------------------------------------------------------------------- ###
    def calc_Ecenter_shift(self):
        """ Eem shift and Eabs shift as a function of temperature """
        
        self.temp:float = np.linspace(prms.tempmin, prms.tempmax, prms.ndiv_temp)
        self.Eabs:float = prms.Eabs0 + ((self.Omegae**2. - self.Omegag**2.)/(self.Omegag**2.)) * prms.kB*self.temp
        self.Eem:float = prms.Eem0 + ((self.Omegag**2. - self.Omegae**2.)/(self.Omegae**2.) +
                                       (8.*(self.Omegag**4.)*self.deltaS)/(self.Omegae**2.*(self.Omegag**2.+self.Omegae**2.)*prms.Eem0)) * prms.kB*self.temp

    ### ----------------------------------------------------------------------------- ###
    def calc_FWHM(self):
        """ Full width half maximum of transitions and its temperature dependence """
        
        self.temp:float = np.linspace(prms.tempmin, prms.tempmax, prms.ndiv_temp)
        self.W0:float = self.Sem * self.Omegag * np.sqrt(8.0*np.log(2.)) / np.sqrt(self.Sabs)
        if ( prms.sw_unit == "nm" ):
            W0min:float = prms.Eem0 - 0.5*self.W0
            W0max:float = prms.Eem0 + 0.5*self.W0
            self.W0 = prms.E2lambda(W0min) - prms.E2lambda(W0max)
        self.W:float = self.W0 * np.sqrt( 1. / np.tanh(self.Omegae/(2.*prms.kB*self.temp)) )
