import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
from scipy import integrate
import sys, os, pathlib, yaml
import cls_subs as subs

prms = subs.subs()

class ph:
    """ Class: Phonon calculations """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor of ph """

        print("* --- Class Phonon --- *")
        if ( prms.sw_phrun ):
            ph.fz_phonon(self)
        if ( not prms.sw_HR == "none" ):
            ph.Huang_Rhys(self)
            ph.plt_specfunc(self)
        print("* --- Finish Class Phonon --- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def fz_phonon(self):
        """ Frozen phonon calculation """
        
        print("* --- Phonon calculation --- *")
        sub.run(["conda activate phonopy"], shell=True)
        sub.run(["phonopy --qe -d --dim='{nd1} {nd2} {nd3}' -c {mat}.scf.in".format(nd1=prms.ndim[0], nd2=prms.ndim[1], nd3 = prms.ndim[2], mat=prms.mat)], shell=True)
        sub.run(["ls -lrt supercell-* > lsout"], shell=True)
        force_command = "phonopy -f "
        data = np.loadtxt("lsout",dtype="str",unpack=True,ndmin=0)
        files = data[len(data)-1]
        for fi in files:
            ncsc = fi.replace("supercell-","").replace(".in","")
            p = pathlib.Path(ncsc)
            if ( not p.is_dir() ):
                sub.run(["mkdir {ncsc}".format(ncsc=ncsc)], shell=True)
                os.chdir(ncsc)
                sub.run(["cat ../header.in ../{fi} >| {mat}-{ncsc}.scf.in".format(fi=fi, mat=prms.mat, ncsc=ncsc)], shell=True)
                """ run scf calculation """
                sub.run(["mpirun -np {nc} {exe}/pw.x < {mat}-{ncsc}.scf.in > {mat}-{ncsc}.scf.out".format(nc=prms.nc, exe=prms.exe, mat=prms.mat, ncsc=ncsc)], shell=True)
                sub.run(["rm work/{mat}.save/wfc*".format(mat=prms.mat)], shell=True)
                os.chdir("../")
            force_command += "{ncsc}/{mat}-{ncsc}.scf.out ".format(mat=prms.mat, ncsc=ncsc)
            
        sub.run([force_command], shell=True)
        pb = pathlib.Path("band.conf")
        if ( p.exists() ):
            sub.run(["phonopy -p band.conf"], shell=True)
        else:
            print("* band.conf does not exist! No plot!")
        print("* --- Finish phonon calculation --- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def dfpt_phonon(self):
        """ Phonon spectrum by DFPT method """

        pass
        
    ### ----------------------------------------------------------------------------- ###
    def Huang_Rhys(self):
        """ obtain Huang_Rhys parameter decomposed by phonon mode """
        
        print("* --- Huang Rhys parameter obtained in Multi-D CCD model --- *")
        print("* use information of Gamma (0,0,0) point")
        p = pathlib.Path("qpoints.yaml")
        if ( not p.exists() ):
            sub.run(["conda activate phonopy"], shell=True)
            sub.run(["phonopy --qpoints='0 0 0' --writedm"], shell=True)
        data = yaml.load(open("qpoints.yaml","r"),Loader=yaml.Loader)
        dynmat_data = data['phonon'][0]['dynamical_matrix']
        for row in dynmat_data:
            vals = np.reshape(row, (-1,2))
            self.dynmat.append(vals[:,0] + vals[:,1] * 1j)
        self.dynmat = np.array(self.dynmat)
        self.eigvals, self.eigvecs = np.linalg.eigh(self.dynmat)
        self.eigvecs = np.reshape(self.eigvecs, (3,len(self.eigvecs[0])//3,len(self.eigvecs)))
        self.freq = np.sqrt(np.abs(self.eigvals.real)) * np.sign(self.eigvals.real)
        self.freq = self.freq * prms.conversion_factor_to_THz * prms.THz2eV

        self.nmode = len(self.freq)
        self.qk = np.zeros(self.nmode)
        self.qk_force = np.zeros(self.nmode)
        (alat, plat, elements, nelems, natm, Rpos_g, volume) = prms.get_POSCAR("POSCAR_"+prms.stateg)
        (alat, plat, elements, nelems, natm, Rpos_e, volume) = prms.get_POSCAR("POSCAR_"+prms.statee)
        mass = []
        for i, nel in enumerate(nelems):
            for j in range(nel):
                mass.append(prms.ELEMS_MASS[elements[i]]*prms.uatm/prms.me)
        if ( prms.sw_HR in {"force","both"} ):
            Force_g = prms.get_FORCE("FORCE_"+prms.stateg)
            Force_e = prms.get_FORCE("FORCE_"+prms.statee)
        
        for i in range(natm):
            for j in range(3):
                if ( prms.sw_HR in {"pos","both"} ):
                    self.qk[:] += np.sqrt(mass[i]) * (Rpos_e[i,j]-Rpos_g[i,j]) * self.eigvecs[j,i,:].real / prms.Bohr
                if ( prms.sw_HR in {"force","both"} ):
                    self.qk_force[:] += (1./((self.freq[:]/prms.Ry)**2.)*np.sqrt(mass[i])) * (Force_e[i,j]-Force_g[i,j]) * self.eigvecs[j,i,:].real
        self.Sk = np.zeros(self.nmode)
        self.Sk_force = np.zeros(self.nmode)
        for i in range(self.nmode):
            self.Sk[i] = self.freq[i] * self.qk[i]**2. / ( 2.*prms.Ry )
            self.Sk_force[i] = self.freq[i] * self.qk_force[i]**2. / ( 2.*prms.Ry )
        energy = np.linspace(prms.emin_ph, prms.emax_ph, prms.ndiv_ph)
        self.Sspec = np.zeros(prms.ndiv_ph)
        self.Sspec_force = np.zeros(prms.ndiv_ph)
        self.Stot = 0.0
        self.Stot_force = 0.0
        for i in range(self.nmode):
            if ( prms.sw_HR in {"pos","both"} ):
                self.Sspec[:] += self.Sk[i] * prms.Gaussian(energy[:]-self.freq[i])
                self.Stot += self.Sk[i]
            if ( prms.sw_HR in {"force","both"} ):
                self.Sspec_force[:] += self.Sk_force[i] * prms.Gaussian(energy[:]-self.freq[i])
                self.Stot_force += self.Sk_force[i]
        print("* Stot from position: ", self.Stot)
        print("* Stot from force: ", self.Stot_force)
        self.DWfac = np.exp(-self.Stot)
        print("* Debye-Waller factor: ", self.DWfac)
        if ( prms.sw_plt_ph ):
            print("* --- Plot S(hbar*omega) from position --- *")
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_xlabel("Phonon energy (meV)")
            ax1.set_ylabel(r"$\mathrm{S(\hbar\omega)}$ (1/meV)")
            ax1.plot(1000.*energy, self.Sspec/1000.,color="black")
            ax2 = ax1.twinx()
            ax2.set_ylabel(r"$\mathrm{S_k}$")
            ax2.bar(1000.*self.freq, self.Sk, color="mediumblue", width=0.15)
            plt.moinorticks_on()
            plt.savefig("Sk_pos.pdf")
            plt.show()
            if ( sw_HR in {"force","both"} ):
                print("* --- Plot S(hbar*omega) from force --- *")
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.set_xlabel("Phonon energy (meV)")
                ax1.set_ylabel(r"$\mathrm{S(\hbar\omega)}$ (1/meV)")
                ax1.plot(1000.*energy, self.Sspec_force/1000.,color="black")
                ax2 = ax1.twinx()
                ax2.set_ylabel(r"$\mathrm{S_k}$")
                ax2.bar(1000.*self.freq, self.Sk_force, color="mediumblue", width=0.15)
                plt.minorticks_on()
                plt.savefig("Sk_force.pdf")
                plt.show()
        print("* --- Finish Huang Rhys parameter --- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def plt_specfunc(self):
        """ plot spectral function """

        print("* --- PLOT SPECTRAL FUNCTION --- *")
        tarr = np.linspace(-prms.tinf, prms.tinf, 2*prms.ndiv_t)
        energy = np.linspace(prms.emin_ph, prms.emax_ph, prms.ndiv_ph)
        integrandS = [self.Sspec * np.exp(-1.j*(energy/prms.hbareV)*t) for t in tarr]
        St = np.array([integrate.simps(integrandS[i,:], energy[:]) for i,t in enumerate(tarr)])
        genfunc = np.array([np.exp(Stt-self.Stot) for Stt in St])
        integrandt = np.array([genfunc*np.exp(1.j*(ene/prms.hbareV)*tarr-prms.gamma_spec*np.abs(tarr)) for ene in enumerate(energy)])
        self.Aspec = (1./(2.*np.pi)) * np.array([integrate.simps(integrandt[i,:],tarr) for i in range(energy)])
        if ( prms.sw_plt_ph ):
            fig = plt.figure()
            plt.xlabel("Energy (eV)")
            plt.ylabel("Intensity (arbitrary unit)")
            plt.plot(energy, self.Aspec, color="mediumblue")
            plt.minortics_on()
            plt.savefig("Spec_MultiD.pdf")
            plt.show()
        else:
            print("* --- NO PLOT --- *")
        print("*")
