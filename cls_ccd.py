import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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
        self.Eabs0_org = prms.Eabs0
        self.Eem0_org = prms.Eem0
        self.EFCg_org = prms.EFCg
        if ( prms.sw_unit in {"eV","nm"} ):
            self.unit = 1.0
        elif ( prms.sw_unit == "cm^-1" ):
            self.unit = prms.eV2cm
        ccd.get_Stokes(self)
        ccd.get_DeltaQ(self)
        if ( prms.sw_eg ):
            ccd.calc_Eeg(self,prms.stateg)
            ccd.calc_Eeg(self,prms.statee)
            ccd.fit_Ecurve(self,"DATA")
        if ( prms.sw_anharm ):
            ccd.fit_Ecurve(self,"CCD")
        if ( not (prms.sw_eg or prms.sw_anharm) ):
            if ( prms.sw_fix == 'curve' ):
                ag = prms.EFCg / (self.dQ_org**2)
                if ( prms.curvature is None ):
                    ae = ag
                else:
                    ae = prms.curvature
                self.deltaQ = (ae*(self.dQ_org**2) + (prms.Eabs0 - (prms.EFCg + prms.Eem0))) / (2.*ae*self.dQ_org)
                self.dQvec = (self.deltaQ/self.dQ_org)*self.dQvec
                self.deltaR = self.deltaQ / self.M
                self.EFCe = ae * self.deltaQ**2
                prms.set_prms_ccd(prms.Eabs0,prms.Eabs0-(ag+ae)*(self.deltaQ**2),ag*(self.deltaQ)**2)
                cfunit = prms.hbar*1.e10*np.sqrt(1./(prms.ep*prms.uatm))
                self.Omegag = cfunit * np.sqrt(2. * prms.EFCg / (self.deltaQ**2.))
                self.Omegae = cfunit * np.sqrt(2. * self.EFCe / (self.deltaQ**2.))
                self.Sem = prms.EFCg / self.Omegag
                self.Sabs = self.EFCe / self.Omegae
                self.deltaS = prms.Eabs0 - prms.Eem0
                self.EZPL = prms.Eem0 + prms.EFCg
        self.temparr = np.linspace(prms.tempmin, prms.tempmax, prms.ndiv_temp)
        ccd.calc_Line_shape(self)
        (self.Eabs, self.Eem) = ccd.calc_Ecenter_shift(self,self.temparr)
        (self.Eabs_temp, self.Eem_temp) = ccd.calc_Ecenter_shift(self,prms.temp)
        (self.W_em, self.W_abs) = ccd.calc_FWHM(self,self.temparr)
        (W_em_temp, W_abs_temp) = ccd.calc_FWHM(self,prms.temp)

        print("*** check the parameters in ccd ***")
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            print("* Absorption energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eabs0*self.unit)
            print("* Emission energy ({sunit}): ".format(sunit=prms.sw_unit), prms.Eem0*self.unit)
            print("* Absorption energy (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), self.Eabs_temp*self.unit)
            print("* Emission energy (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), self.Eem_temp*self.unit)
            print("* EFCg = Ee - Eg ({sunit}): ".format(sunit=prms.sw_unit), prms.EFCg*self.unit)
            print("* EFCe = Eg* - Ee* ({sunit}): ".format(sunit=prms.sw_unit), self.EFCe*self.unit)
            print("* Stokes shift (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.deltaS*self.unit)
            print("* Zero phonon photoemission line ({sunit}): ".format(sunit=prms.sw_unit), self.EZPL*self.unit)
            print("* FWHM_em (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.W0_em*self.unit)
            print("* FWHM_abs (@ 0K) ({sunit}): ".format(sunit=prms.sw_unit), self.W0_abs*self.unit)
            print("* FWHM_em (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), W_em_temp*self.unit)
            print("* FWHM_abs (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), W_abs_temp*self.unit)
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3] (sqrt(amu)*ang): ", self.dQvec)
            print("* hbar*Omegag ({sunit}): ".format(sunit=prms.sw_unit), self.Omegag*self.unit)
            print("* hbar*Omegae ({sunit}): ".format(sunit=prms.sw_unit), self.Omegae*self.unit)
        elif ( prms.sw_unit == "nm" ):
            print("* Absorption energy (nm): ", prms.E2lambda(prms.Eabs0))
            print("* Emission energy (nm): ", prms.E2lambda(prms.Eem0))
            print("* Absorption energy (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), prms.E2lambda(self.Eabs_temp))
            print("* Emission energy (@ {temp}K) ({sunit}): ".format(temp=prms.temp,sunit=prms.sw_unit), prms.E2lambda(self.Eem_temp))
            print("* EFCg = Ee - Eg (nm): ", prms.E2lambda(prms.EFCg))
            print("* EFCe = Eg* - Ee* (nm): ", prms.E2lambda(self.EFCe))
            print("* Stokes shift (@ 0K) (nm): ", prms.E2lambda(self.deltaS))
            print("* Zero phonon photoemission line (nm): ", prms.E2lambda(self.EZPL))
            print("* FWHM_em (@ 0K) (nm): ", self.W0_em)
            print("* FWHM_abs (@ 0K) (nm): ", self.W0_abs)
            print("* FWHM_em (@ {temp}K) (nm): ".format(temp=prms.temp), W_em_temp)
            print("* FWHM_abs (@ {temp}K) (nm): ".format(temp=prms.temp), W_abs_temp)
            print("* DeltaR (ang): ", self.deltaR)
            print("* DeltaQ (sqrt(amu)*ang): ", self.deltaQ)
            print("* DeltaQ[1:3]: (sqrt(amu)*ang)", self.dQvec)
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

        print("*** PLOT 1D-CCD ***")
        plt.figure(figsize=(6,6.5))
        ccd.plt_1DCCD(self)
        #plt.minorticks_on()
        #plt.savefig("1DCCD_{state}.pdf".format(state=prms.statee))
        #plt.show()
        if ( prms.sw_eg or prms.sw_anharm ):
            ccd.plt_1DCCD(self,sw_anharm=True)
        plt.minorticks_on()
        plt.savefig("1DCCD_{state}.pdf".format(state=prms.statee))
        plt.show()
        
        print("*** PLOT LINE SHAPE ***")
        size = 30
        plt.ylabel("Intensity (arbitrary unit)")
        plt.axhline(0, lw=0.5, c="black", linestyle="dashed")
        if ( not prms.include == None ):
            energy_exp,int_exp = np.loadtxt(prms.include,dtype='float',unpack=True,ndmin=0)
            energy_exp = prms.fac_include * energy_exp
            plt.scatter(energy_exp,int_exp,c="white",edgecolor="mediumblue",marker="o",s=80)
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            if ( prms.sw_spec_temp ):
                plt.xlabel("Energy ({sunit})".format(sunit=prms.sw_unit))
                plt.xlim(self.unit*prms.emin_ccd,self.unit*prms.emax_ccd)
                plt.plot(self.unit*(self.energy-prms.Eabs0+self.Eabs_temp), self.Labs, lw=1.0, c="black")
                plt.scatter(self.Eabs_temp*self.unit,0.0,c="black",s=size,marker="o")
                plt.plot(self.unit*(self.energy-prms.Eem0+self.Eem_temp), self.Lem, lw=1.0, c="red")
                plt.scatter(self.Eem_temp*self.unit,0.0,c="red",s=size,marker="o")
            else:
                plt.xlabel("Energy ({sunit})".format(sunit=prms.sw_unit))
                plt.plot(self.unit*self.energy, self.Labs, lw=1.0, c="black")
                plt.scatter(prms.Eabs0*self.unit,0.0,c="black",s=size,marker="o")
                plt.plot(self.unit*self.energy, self.Lem, lw=1.0, c="red")
                plt.scatter(prms.Eem0*self.unit,0.0,c="red",s=size,marker="o")
        elif ( prms.sw_unit == "nm" ):
            if ( prms.sw_spec_temp ):
                plt.xlabel("Wave length (nm)")
                plt.plot(prms.E2lambda(self.energy-prms.Eabs0+self.Eabs_temp), 0.9*self.Labs, lw=1.0, c="black")
                plt.plot(prms.E2lambda(self.energy-prms.Eem0+self.Eem_temp), 0.9*self.Lem, lw=1.0, c="red")
            else:
                plt.xlabel("Wave length (nm)")
                plt.plot(prms.E2lambda(self.energy), self.Labs, lw=1.0, c="black")
                plt.plot(prms.E2lambda(self.energy), self.Lem, lw=1.0, c="red")
        plt.minorticks_on()
        plt.savefig("Spectrum.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF E_EM & E_ABS ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Emission energy ({sunit})".format(sunit=prms.sw_unit))
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            plt.plot(self.temparr, self.unit*self.Eabs, lw=1.0, c="black")
            plt.plot(self.temparr, self.unit*self.Eem, lw=1.0, c="red")
        elif ( prms.sw_unit == "nm" ):
            plt.plot(self.temparr, prms.E2lambda(self.Eabs), lw=1.0, c="black")
            plt.plot(self.temparr, prms.E2lambda(self.Eem), lw=1.0, c="red")
        plt.minorticks_on()
        plt.savefig("Epeak_Temp.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF STOKES SHIFT ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Stokes shift ({sunit})".format(sunit=prms.sw_unit))
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            plt.plot(self.temparr, self.unit*(self.Eabs-self.Eem), lw=1.0, c="black")
        elif ( prms.sw_unit == "nm" ):
            plt.plot(self.temparr, prms.E2lambda(self.Eabs)-prms.E2lambda(self.Eem), lw=1.0, c="black")
        plt.minorticks_on()
        plt.savefig("Stokes_Temp.pdf")
        plt.show()
                
        print("*** PLOT TEMPERATURE DEPENDENCE OF FWHM ***")
        plt.xlabel("Temperature (K)")
        plt.ylabel("FWHM ({sunit})".format(sunit=prms.sw_unit))
        plt.plot(self.temparr, self.unit*self.W_abs, lw=1.0, c="black")
        plt.plot(self.temparr, self.unit*self.W_em, lw=1.0, c="red")
        plt.minorticks_on()
        plt.savefig("FWHM_Temp.pdf")
        plt.show()

    ### ----------------------------------------------------------------------------- ###
    def plt_1DCCD(self, sw_anharm = False):
        """ plot 1D-CCD """

        def Ecurve(x,a,Q0,E0):
            return a*(x-Q0)**2 + E0

        if ( not sw_anharm ):
            if ( prms.sw_fix == 'vertex' ):
                Eg0 = 0.0
                # Ee0 = prms.EFCg + prms.Eem0
                Ee0 = self.Eem_2dim+self.EFCg_2dim
                ag = self.EFCg_2dim / (self.dQ_org**2)
                ae = self.EFCe_2dim / (self.dQ_org**2)
                Qg0 = 0.0
                Qe0 = self.dQ_org
            elif ( prms.sw_fix == 'curve' ):
                ag = self.EFCg_2dim / (self.deltaQ**2)
                if ( prms.curvature is None ):
                    ae = ag
                else:
                    ae = prms.curvature
                Qg0 = 0.0
                Qe0 = self.deltaQ
                Eg0 = 0.0
                Ee0 = self.Eem_2dim+self.EFCg_2dim
            else:
                print("*** ERROR in plt_1DCCD: sw_fix should be 'vertex' or 'curve'!!!")
                sys.exit()
                
            plt.xlim(-0.5*self.dQ_org,1.5*self.dQ_org)
            plt.ylim(-0.1, prms.emax_ccd)
            plt.xlabel(r"$\Delta Q$ (amu$^{1/2}\cdot\mathrm{\AA}$)")
            plt.ylabel("Energy (eV)")
            xQ = np.linspace(-0.5*self.dQ_org,1.5*self.dQ_org,1000)
            yEg = Ecurve(xQ,ag,Qg0,Eg0)
            yEe = Ecurve(xQ,ae,Qe0,Ee0)
            Qarr = [0.0, 0.0, self.dQ_org, self.dQ_org]
            Earr = [0.0, self.Eabs0_org, self.EFCg_org, self.EFCg_org+self.Eem0_org]
            plt.scatter(Qarr, Earr, color="white", marker="o", edgecolor="mediumblue", s=80)
            plt.plot(xQ,yEg,color="mediumblue",lw=1.0,linestyle='dashed')
            plt.plot(xQ,yEe,color="mediumblue",lw=1.0,linestyle='dashed')
            plt.plot([0.0,0.0],[-0.1,self.Eabs0_org],color="black",linestyle="dotted",lw=0.5)
            plt.plot([Qe0,Qe0],[-0.1,Ee0],color="black",linestyle="dotted",lw=0.5)
            plt.plot([-0.5*self.dQ_org,self.deltaQ],[Ee0,Ee0],color="black",linestyle="dotted",lw=0.5)
            plt.plot([-0.5*self.dQ_org,self.deltaQ],[Ecurve(Qe0,ag,Qg0,Eg0),Ecurve(Qe0,ag,Qg0,Eg0)],color="black",linestyle="dotted",lw=0.5)
            plt.plot([-0.5*self.dQ_org,1.5*self.dQ_org],[0.0,0.0],color="black",linestyle="dotted",lw=0.5)
            plt.plot([-0.5*self.dQ_org,0.0],[self.Eabs0_org,self.Eabs0_org],color="black",linestyle="dotted",lw=0.5)

        else:
            if ( prms.sw_anharm ):
                fn_g = "CCD_{state}".format(state=prms.stateg)
                fn_e = "CCD_{state}".format(state=prms.statee)                
            if ( prms.sw_eg ):
                fn_g = "DATA_{state}".format(state=prms.stateg)
                fn_e = "DATA_{state}".format(state=prms.statee)

            plt.xlim(-0.5*self.dQ_org,1.5*self.dQ_org)
            plt.ylim(-0.1, prms.emax_ccd)
            plt.xlabel(r"$\Delta Q$ (amu$^{1/2}\cdot\mathrm{\AA}$)")
            plt.ylabel("Energy (eV)")
            Q, Etot_g = np.loadtxt(fn_g,dtype='float',unpack=True,ndmin=0)
            x = np.linspace(-0.5*self.dQ_org+Q[0],0.5*self.dQ_org+Q[len(Q)-1],1000)
            weight = np.ones(len(Q))
            coef_egg = np.polyfit(Q, Etot_g, 3, w=weight)
            Efit_egg = np.poly1d(coef_egg)(x)
            Qmin_egg = (-coef_egg[1]+np.sqrt(coef_egg[1]**2.-3.*coef_egg[0]*coef_egg[2]))/(3.*coef_egg[0])
            Q, Etot_e = np.loadtxt(fn_e,dtype='float',unpack=True,ndmin=0)
            coef_ege = np.polyfit(Q, Etot_e, 3, w=weight)
            Efit_ege = np.poly1d(coef_ege)(x)
            Qmin_ege = (-coef_ege[1]+np.sqrt(coef_ege[1]**2.-3.*coef_ege[0]*coef_ege[2]))/(3.*coef_ege[0])
            plt.scatter(Q,(Etot_g-min(Etot_g))*prms.Ry,marker="o",label="data",edgecolor="crimson",color="white",s=50)
            plt.plot(x,(Efit_egg-min(Etot_g))*prms.Ry,linestyle="dashed",color="red")
            plt.scatter(Qmin_egg, (np.poly1d(coef_egg)(Qmin_egg)-min(Etot_g))*prms.Ry, marker="*", color="darkblue", s=60)
            plt.scatter(Q,(Etot_e-min(Etot_g))*prms.Ry,marker="o",label="data",edgecolor="crimson",color="white",s=50)
            plt.plot(x,(Efit_ege-min(Etot_g))*prms.Ry,linestyle="dashed",color="red")
            plt.scatter(Qmin_ege, (np.poly1d(coef_ege)(Qmin_ege)-min(Etot_g))*prms.Ry, marker="*", color="darkblue", s=60)
            
    ### ----------------------------------------------------------------------------- ###
    def calc_Eeg(self, state):
        """ calculate total enegies of intermediate states """
        
        print("* --- gernerate POSCAR in between ground & excited states --- *")
        (alat_g,plat_g,elements_g,nelems_g,natm_g,pos_g,volume) = prms.get_POSCAR("POSCAR_"+prms.stateg)
        (alat_e,plat_e,elements_e,nelems_e,natm_e,pos_e,volume) = prms.get_POSCAR("POSCAR_"+prms.statee)
        Qarr = np.linspace(-prms.dQ, 1.+prms.dQ, prms.ndiv_eg)
        for i, Q in enumerate(Qarr):
            Rmid = np.zeros((len(pos_g),3))
            Rmid[:,:] = Q * (pos_e[:,:] - pos_g[:,:]) + pos_g[:,:]
            fn_eg = "POSCAR_eg{ic}".format(ic=i+1)
            sub.run(["head -8 POSCAR_{state} > {fn}".format(state=prms.stateg, fn=fn_eg)], shell=True)
            strcoord = ""
            for j in range(len(Rmid)):
                strcoord += " {posx}  {posy}  {posz} \n".format(posx=Rmid[j,0], posy=Rmid[j,1], posz=Rmid[j,2])
            with open(fn_eg, "a") as f:
                f.write(strcoord)
            print("* Q={Q}*deltaQ {fn} created!".format(Q=Q, fn=fn_eg))
        print("* Finish!")
        print("*")

        fn = "DATA_{state}".format(state=state)
        fp = pathlib.Path(fn)
        if ( not fp.exists() ):
            string = "### 1:Q((amu)^(1/2)*ang)   2:Etot(Ry) \n"
            with open(fn,"w") as f:
                f.write(string)
        for i, qc in enumerate(Qarr):
            ndir = "eg{ic}".format(ic=i+1)
            pd = pathlib.Path("{ndir}/{state}".format(ndir=ndir, state=state))
            if ( not pd.is_dir() ):
                sub.run(["mkdir -p {ndir}/{state}".format(ndir=ndir, state=state)], shell=True)
                os.chdir("{ndir}/{state}".format(ndir=ndir, state=state))
                sub.run(["cp ../../header_{state}.in {mat}-eg{ic}.scf.in".format(state=state, mat=prms.mat, ic=i+1)], shell=True)
                (alat, plat, elements, nelems, natm, pos, volume) = prms.get_POSCAR("../../POSCAR_eg{ic}".format(ic=i+1))
                ic = 0
                string = ""
                for j, nele in enumerate(nelems):
                    for k in range(nele):
                        string += " {0:}  {1:.10f}  {2:.10f}  {3:.10f}  \n".format(elements[j],pos[ic,0],pos[ic,1],pos[ic,2])
                        ic += 1
                with open("{mat}-eg{ic}.scf.in".format(mat=prms.mat,ic=i+1),"a") as f:
                    f.write(string)
                """ run scf calculation """
                sub.run(["mpirun -np {nc} {exe}/pw.x < {mat}-eg{ic}.scf.in > {mat}-eg{ic}.scf.out".format(nc=prms.nc, exe=prms.exe, ic=i+1, mat=prms.mat)], shell=True)
                sub.run(["rm work/{mat}.save/wfc*".format(mat=prms.mat)], shell=True)
                sub.run(["grep ! {mat}-eg{ic}.scf.out > grep.out".format(mat=prms.mat, ic=i+1)], shell=True)
                data = np.loadtxt("grep.out",dtype="str",unpack=True,ndmin=0)
                os.chdir("../../")
                string = "{Q:.10f}   {dat:.10f} \n".format(Q=qc*self.deltaQ,dat=float(data[4]))
                with open(fn, "a") as f:
                    f.write(string)

    ### ----------------------------------------------------------------------------- ###
    def fit_Ecurve(self, fn):
        """ fiting Energy curve by polynomial """

        print("* --- Fitting energy curve --- *")
        fng = fn+"_"+prms.stateg
        fne = fn+"_"+prms.statee
        Qg, Etotg = np.loadtxt(fng,dtype='float',unpack=True,ndmin=0)
        Qe, Etote = np.loadtxt(fne,dtype='float',unpack=True,ndmin=0)
        xg = np.linspace(-0.1+Qg[0],0.1+Qg[len(Qg)-1],1000)
        xe = np.linspace(-0.1+Qe[0],0.1+Qe[len(Qe)-1],1000)
        weightg = np.ones(len(Qg))
        weightg_add = np.append(weightg,1.e3)
        weighte = np.ones(len(Qe))
        Qg_add = np.append(Qg,0.)
        Etotg_add = np.append(Etotg,min(Etotg))
        coefg_3dim = np.polyfit(Qg_add, (Etotg_add-min(Etotg_add))*prms.Ry, 3, w=weightg_add)
        Efitg_3dim = np.poly1d(coefg_3dim)(xg)
        coefg_2dim = np.polyfit(Qg_add, (Etotg_add-min(Etotg_add))*prms.Ry, 2, w=weightg_add)
        Efitg_2dim = np.poly1d(coefg_2dim)(xg)
        coefe_3dim = np.polyfit(Qe, (Etote-min(Etote))*prms.Ry, 3, w=weighte)
        Efite_3dim = np.poly1d(coefe_3dim)(xe)
        coefe_2dim = np.polyfit(Qe, (Etote-min(Etote))*prms.Ry, 2, w=weighte)
        Efite_2dim = np.poly1d(coefe_2dim)(xe)

        """ consider Anharmonic effects """
        cfunit = prms.hbar*1.e10*np.sqrt(1./(prms.ep*prms.uatm))
        Qming_3dim = (-coefg_3dim[1]+np.sqrt(coefg_3dim[1]**2.-3.*coefg_3dim[0]*coefg_3dim[2]))/(3.*coefg_3dim[0])
        Qmine_3dim = (-coefe_3dim[1]+np.sqrt(coefe_3dim[1]**2.-3.*coefe_3dim[0]*coefe_3dim[2]))/(3.*coefe_3dim[0])
        self.Omegag_3dim = cfunit*np.sqrt(6.*coefg_3dim[0]*Qming_3dim+2.*coefg_3dim[1])
        self.EFCg_3dim = np.poly1d(coefg_3dim)(Qmine_3dim)-np.poly1d(coefg_3dim)(0.)
        self.Eabs_3dim = np.poly1d(coefe_3dim)(0.)+(min(Etote)-min(Etotg))*prms.Ry
        self.Eem_3dim = np.poly1d(coefe_3dim)(Qmine_3dim)+(min(Etote)-min(Etotg))*prms.Ry-self.EFCg_3dim
        self.Sem_3dim = self.EFCg_3dim / self.Omegag_3dim
        self.Omegae_3dim = cfunit*np.sqrt(6.*coefe_3dim[0]*Qmine_3dim+2.*coefe_3dim[1])
        self.EFCe_3dim = np.poly1d(coefe_3dim)(0.)-np.poly1d(coefe_3dim)(Qmine_3dim)
        self.Sabs_3dim = self.EFCe_3dim / self.Omegae_3dim
        self.dQvec_3dim = ( Qmine_3dim / self.dQ_org) * self.dQvec
        self.deltaQ_3dim = Qmine_3dim
        self.deltaR_3dim = self.deltaQ_3dim / self.M
        self.deltaS_3dim = self.Eabs_3dim - ( self.Eem_3dim - self.EFCg_3dim )
        self.EZPL_3dim = self.Eem_3dim
        print("* Fitted by cubic equation: ")
        print("* ground state *")
        print(coefg_3dim)
        print("* hbar*Omegag (eV) [cubic]: sqrt({cube}Q+{square})".format(cube=(cfunit**2.)*6.*coefg_3dim[0], square=(cfunit**2.)*2.*coefg_3dim[1]))
        print("* Qming (amu^1/2 ang): {Qmin}".format(Qmin=Qming_3dim))
        print("* hbar*Omegag (eV) @ Q=Qming: {cubic}".format(cubic=self.Omegag_3dim))
        print("*")
        print("* excited state *")
        print(coefe_3dim)
        print("* hbar*Omegae (eV) [cubic]: sqrt({cube}Q+{square})".format(cube=(cfunit**2.)*6.*coefe_3dim[0], square=(cfunit**2.)*2.*coefe_3dim[1]))
        print("* Qmine (amu^1/2 ang): {Qmin}".format(Qmin=Qmine_3dim))
        print("* hbar*Omegae (eV) @ Q=Qmine: {cubic}".format(cubic=self.Omegae_3dim))
        print("*")
    
        """ Harmonic approximation """
        self.Omegag_2dim = cfunit*np.sqrt(2.*coefg_2dim[0])
        Qmine_2dim = - coefe_2dim[1] / (2.*coefe_2dim[0])
        Qming_2dim = - coefg_2dim[1] / (2.*coefg_2dim[0])
        self.EFCg_2dim = np.poly1d(coefg_2dim)(Qmine_2dim)-np.poly1d(coefg_2dim)(0.)
        self.Eabs_2dim = np.poly1d(coefe_2dim)(0.)+(min(Etote)-min(Etotg))*prms.Ry
        self.Eem_2dim = np.poly1d(coefe_2dim)(Qmine_2dim)+(min(Etote)-min(Etotg))*prms.Ry-self.EFCg_2dim
        self.Sem_2dim = self.EFCg_2dim / self.Omegag_2dim
        self.Omegae_2dim = cfunit*np.sqrt(2.*coefe_2dim[0])
        self.EFCe_2dim = np.poly1d(coefe_2dim)(0.)-np.poly1d(coefe_2dim)(Qmine_2dim)
        self.Sabs_2dim = self.EFCe_2dim / self.Omegae_2dim
        self.dQvec_2dim = ( Qmine_2dim / self.dQ_org) * self.dQvec
        self.deltaQ_2dim = Qmine_2dim 
        self.deltaR_2dim = self.deltaQ_2dim / self.M
        self.deltaS_2dim = self.Eabs_2dim - ( self.Eem_2dim - self.EFCg_2dim )
        self.EZPL_2dim = self.Eem_2dim 
        print("* Fitted by quadratic equation: ")
        print("* ground state *")
        print(coefg_2dim)
        print("* hbar*Omegag (eV) [quadratic]: {square}".format(square=self.Omegag_2dim))
        print("* Qming (amu^1/2 ang): {Qmin}".format(Qmin=Qming_2dim))
        print("* excited state *")
        print(coefe_2dim)
        print("* hbar*Omegae (eV) [quadratic]: {square}".format(square=self.Omegae_2dim))
        print("* Qmine (amu^1/2 ang): {Qmin}".format(Qmin=Qmine_2dim))

        if ( prms.dim_fit == 2 ):
            self.Omegag = self.Omegag_2dim
            self.Omegae = self.Omegae_2dim
            prms.set_prms_ccd(self.Eabs_2dim,self.Eem_2dim-self.EFCg_2dim,self.EFCg_2dim)
            self.EFCe = self.EFCe_2dim
            self.Sabs = self.Sabs_2dim
            self.Sem = self.Sem_2dim
            self.dQvec = self.dQvec_2dim
            self.deltaQ = self.deltaQ_2dim
            self.deltaR = self.deltaR_2dim
            self.deltaS = self.deltaS_2dim
            self.EZPL = self.EZPL_2dim
        else:
            self.Omegag = self.Omegag_3dim
            self.Omegae = self.Omegae_3dim
            prms.set_prms_ccd(self.Eabs_3dim,self.Eem_3dim-self.EFCg_3dim,self.EFCg_3dim)
            self.EFCe = self.EFCe_3dim
            self.Sabs = self.Sabs_3dim
            self.Sem = self.Sem_3dim
            self.dQvec = self.dQvec_3dim
            self.deltaQ = self.deltaQ_3dim
            self.deltaR = self.deltaR_3dim
            self.deltaS = self.deltaS_3dim
            self.EZPL = self.EZPL_3dim
        
        if prms.sw_plt_ccd:
            plt.xlabel(r"Q ($\mathrm{amu}^{1/2} \cdot \AA$)")
            plt.ylabel("Energy (eV)")
            # plt.xlim(-0.05*self.dQ_org, 1.05*self.dQ_org)
            plt.ylim(-0.1,prms.emax_ccd)
            size = 50
            """ ground state """
            plt.scatter(Qg,(Etotg-min(Etotg))*prms.Ry,marker="o",edgecolor="darkred",color="white",s=size)
            if prms.sw_2dim:
                plt.plot(xg,Efitg_2dim,linestyle="dashed",color="mediumblue")
                plt.scatter(Qming_2dim, np.poly1d(coefg_2dim)(Qming_2dim), marker="*", color="darkblue", s=2*size)
            plt.plot(xg,Efitg_3dim,linestyle="dashed",color="coral")
            plt.scatter(Qming_3dim, np.poly1d(coefg_3dim)(Qming_3dim), marker="*", color="crimson", s=2*size)
            """ excited state """
            plt.scatter(Qe,(Etote-min(Etotg))*prms.Ry,marker="o",label="data",edgecolor="darkred",color="white",s=size)
            if prms.sw_2dim:
                pass
                plt.plot(xe,Efite_2dim+(min(Etote)-min(Etotg))*prms.Ry,linestyle="dashed",color="mediumblue",label="fit 2dim")
                plt.scatter(Qmine_2dim, np.poly1d(coefe_2dim)(Qmine_2dim)+min(Etote)*prms.Ry, marker="*", color="darkblue", s=2*size)
            plt.plot(xe,Efite_3dim+(min(Etote)-min(Etotg))*prms.Ry,linestyle="dashed",color="coral",label="fit 3dim")
            plt.scatter(Qmine_3dim, np.poly1d(coefe_3dim)(Qmine_3dim)+(min(Etote)-min(Etotg))*prms.Ry, marker="*", color="crimson", s=2*size)
            plt.legend(fontsize=12)
            plt.minorticks_on()
            plt.savefig("Ecurve_fit_{fn}_{state}.pdf".format(fn=fn,state=prms.statee))
            plt.show()
        
    ### ----------------------------------------------------------------------------- ###
    def get_Stokes(self):
        """ Stokes shift, Franck-Condon parameter, and zero-phonon energy """
        
        self.deltaS = prms.Eabs0 - prms.Eem0
        self.EFCe = self.deltaS - prms.EFCg
        self.EZPL = prms.Eem0 + prms.EFCg

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
            
        mass = [prms.ELEMS_MASS[ele]*prms.uatm/prms.me for ele in elements_g]
        pos_g = alat_g * np.dot( pos_g, plat_g )
        pos_e = alat_e * np.dot( pos_e, plat_e )
        deltaQ_sq = 0.0
        deltaR_sq = 0.0
        dQvec_sq = np.zeros(3)
        count:int = 0
        for i, nel in enumerate(nelems_g):
            for j in range(nel):
                deltaQ_sq += mass[i] * ( sum([(pos_e[count,k]-pos_g[count,k])**2. for k in range(3)]) )
                deltaR_sq += sum([(pos_e[count,k]-pos_g[count,k])**2. for k in range(3)])
                dQvec_sq[:] += mass[i] * (pos_e[count,:]-pos_g[count,:])**2.
                count += 1
        self.deltaQ = np.sqrt(prms.me/prms.uatm) * prms.Bohr * np.sqrt(deltaQ_sq)
        self.dQ_org = self.deltaQ
        self.dQvec = np.sqrt(prms.me/prms.uatm) * prms.Bohr * np.array([np.sqrt(dQvec_sq[i]) for i in range(3)])
        self.deltaR = prms.Bohr * np.sqrt(deltaR_sq)
        self.M = deltaQ_sq / deltaR_sq

        cfunit = prms.hbar*1.e10*np.sqrt(1./(prms.ep*prms.uatm))
        self.Omegag = cfunit*np.sqrt(2. * prms.EFCg / (self.deltaQ**2.)) 
        self.Omegae = cfunit*np.sqrt(2. * self.EFCe / (self.deltaQ**2.)) 
        self.Sabs = self.EFCe / self.Omegae
        self.Sem = prms.EFCg / self.Omegag
        
    ### ----------------------------------------------------------------------------- ###
    def calc_Line_shape(self):
        """ Line shape of optical spectra """

        room = 0.3
        if ( prms.sw_unit in {"eV","cm^-1"} ):
            self.energy = np.linspace(prms.emin_ccd-room, prms.emax_ccd+room, prms.ndiv_e)
        else:
            self.energy = np.linspace(prms.E2lambda(prms.emax_ccd+room),prms.E2lambda(prms.emin_ccd-room),prms.ndiv_e)
            self.energy = prms.lambda2E(self.energy)
        self.Lem = np.zeros(prms.ndiv_e)
        self.Labs = np.zeros(prms.ndiv_e)
        for n in range(prms.nmax):
            self.Lem[:] += ( np.exp(-self.Sem)*self.Sem**float(n) / float(np.math.factorial(n)) ) * prms.Lorentzian(self.EZPL - float(n)*self.Omegag - self.energy[:])
            self.Labs[:] += ( np.exp(-self.Sabs)*self.Sabs**float(n) / float(np.math.factorial(n)) ) * prms.Lorentzian(self.EZPL + float(n)*self.Omegae - self.energy[:])
        self.Lem = prms.I0 * self.Lem / max(self.Lem)
        self.Labs = prms.I0 * self.Labs / max(self.Labs)

    ### ----------------------------------------------------------------------------- ###
    def calc_Ecenter_shift(self, temperature):
        """ Eem shift and Eabs shift as a function of temperature """
        
        Eabs_temp = prms.Eabs0 + ((self.Omegae**2. - self.Omegag**2.)/(self.Omegag**2.)) * prms.kB*temperature
        Eem_temp = prms.Eem0 + ((self.Omegag**2. - self.Omegae**2.)/(self.Omegae**2.) +
                                      (8.*(self.Omegag**4.)*self.deltaS)/(self.Omegae**2.*(self.Omegag**2.+self.Omegae**2.)*prms.Eem0)) * prms.kB*temperature
        return ( Eabs_temp, Eem_temp )
        
    ### ----------------------------------------------------------------------------- ###
    def calc_FWHM(self, temperature):
        """ Full width half maximum of transitions as a function of temperature """
        
        self.W0_em = self.Sem * self.Omegag * np.sqrt(8.0*np.log(2.)) / np.sqrt(self.Sabs)
        self.W0_abs = self.Sabs * self.Omegae * np.sqrt(8.0*np.log(2.)) / np.sqrt(self.Sem)
        if ( prms.sw_unit == "nm" ):
            W0minem = prms.Eem0 - 0.5*self.W0_em
            W0maxem = prms.Eem0 + 0.5*self.W0_em
            self.W0_em = prms.E2lambda(W0minem) - prms.E2lambda(W0maxem)
            W0minabs = prms.Eabs0 - 0.5*self.W0_abs
            W0maxabs = prms.Eabs0 + 0.5*self.W0_abs
            self.W0_abs = prms.E2lambda(W0minabs) - prms.E2lambda(W0maxabs)
        W_em_temp = self.W0_em * np.sqrt( 1. / np.tanh(self.Omegae/(2.*prms.kB*temperature)) )
        W_abs_temp = self.W0_abs * np.sqrt( 1. / np.tanh(self.Omegag/(2.*prms.kB*temperature)) )
        return ( W_em_temp, W_abs_temp )
