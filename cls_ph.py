import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
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
            ph.phonon(self)
        if ( not prms.sw_HR == "none" ):
            ph.Huang_Rhys(self)
        print("* --- Finish Class Phonon --- *")
        print("*")

    ### ----------------------------------------------------------------------------- ###
    def phonon(self):
        """ Frozen phonon calculation """
        
        print("* --- Phonon calculation --- *")
        sub.run(["conda activate phonopy"], shell=True)
        sub.run(["phonopy --qe -d --dim='{nd1} {nd2} {nd3}' -c {mat}.scf.in".format(nd1=prms.ndim[0], nd2=prms.ndim[1], nd3 = prms.ndim[2], mat=prms.mat)], shell=True)
        sub.run(["ls -lrt supercell-* > lsout"], shell=True)
        force_command:str = "phonopy -f "
        data:str = np.loadtxt("lsout",dtype="str",unpack=True,ndmin=0)
        files:str = data[len(data)-1]
        for fi in files:
            ncsc:str = fi.replace("supercell-","").replace(".in","")
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
        self.freq:float = np.sqrt(np.abs(self.eigvals.real)) * np.sign(self.eigvals.real)
        self.freq = self.freq * prms.conversion_factor_to_THz * prms.THz2eV

        self.nmode:int = len(self.freq)
        self.qk:float = np.zeros(self.nmode)
        self.qk_force:float = np.zeros(self.nmode)
        (alat, plat, elements, nelems, natm, Rpos_g, volume) = prms.get_POSCAR("POSCAR_"+prms.stateg)
        (alat, plat, elements, nelems, natm, Rpos_e, volume) = prms.get_POSCAR("POSCAR_"+prms.statee)
        mass:float = []
        for i, nel in enumerate(nelems):
            for j in range(nel):
                mass.append(prms.ELEMS_MASS[elements[i]]*prms.uatm/prms.me)
        if ( sw_qk in {"force","both"} ):
            Force_g:float = prms.get_FORCE("FORCE_"+prms.stateg)
            Force_e:float = prms.get_FORCE("FORCE_"+prms.statee)
        
        for i in range(natm):
            for j in range(3):
                if ( prms.sw_qk in {"pos","both"} ):
                    self.qk[:] += np.sqrt(mass[i]) * (Rpos_e[i,j]-Rpos_g[i,j]) * self.eigvecs[j,i,:].real / prms.Bohr
                if ( prms.sw_qk in {"force","both"} ):
                    self.qk_force[:] += (1./((self.freq[:]/prms.Ry)**2.)*np.sqrt(mass[i])) * (Force_e[i,j]-Force_g[i,j]) * self.eigvecs[j,i,:].real
        self.Sk:float = np.zeros(self.nmode)
        self.Sk_force:float = np.zeros(self.nmode)
        for i in range(self.nmode):
            self.Sk[i]:float = self.freq[i] * self.qk[i]**2. / ( 2.*prms.Ry )
            self.Sk_force[i]:float = self.freq[i] * self.qk_force[i]**2. / ( 2.*prms.Ry )
        energy:float = np.linspace(prms.emin_ph, prms.emax_ph, prms.ndiv_ph)
        self.Sspec:float = np.zeros(prms.ndiv_ph)
        self.Sspec_force:float = np.zeros(prms.ndiv_ph)
        self.Stot:float = 0.0
        self.Stot_force:float = 0.0
        for i in range(self.nmode):
            if ( prms.sw_qk in {"pos","both"} ):
                self.Sspec[:] += self.Sk[i] * prms.Gaussian(energy[:]-self.freq[i])
                self.Stot += self.Sk[i]
            if ( prms.sw_qk in {"force","both"} ):
                self.Sspec_force[:] += self.Sk_force[i] * prms.Gaussian(energy[:]-self.freq[i])
                self.Stot_force += self.Sk_force[i]
        print("* Stot from position: ", self.Stot)
        print("* Stot from force: ", self.Stot_force)
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
            plt.show()
            if ( sw_qk in {"force","both"} ):
                print("* --- Plot S(hbar*omega) from force --- *")
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.set_xlabel("Phonon energy (meV)")
                ax1.set_ylabel(r"$\mathrm{S(\hbar\omega)}$ (1/meV)")
                ax1.plot(1000.*energy, self.Sspec_force/1000.,color="black")
                ax2 = ax1.twinx()
                ax2.set_ylabel(r"$\mathrm{S_k}$")
                ax2.bar(1000.*self.freq, self.Sk_force, color="mediumblue", width=0.15)
                plt.show()
        print("* --- Finish Huang Rhys parameter --- *")
        print("*")
