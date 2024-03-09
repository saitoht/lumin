<<<<<<< HEAD
### 9th, Mar., 2024   H. Saito
=======
### 12th, Feb., 2024   H. Saito
>>>>>>> 2482b5ac0daf49ca16f8c4b31e4e33b5393f7357
### See "A. Alkauskas, B. B. Buckley, D. D. Awschalom, & C. G. Van de Walle,
###      First-principels theory of the luminescence lineshape for the triplet transition in diamond NV centres, New J. Phys. 16 (2014) 073026."
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import sys, os, pathlib, yaml
import cls_subs as subs

class ph:
    """ Class: Phonon calculations """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        const = subs.subs()
        prms = subs.get_prms()

    ### ----------------------------------------------------------------------------- ###
    def run_phonon():
        """ execute phonon calculations """    
        if ( prms.sw_phrun ):
            phonon()
        if ( not prms.sw_HR == "none" ):
            Huang_Rhys()

    ### ----------------------------------------------------------------------------- ###
    def phonon(self):
        """ Frozen phonon calculation """
        print("* --- Phonon calculation --- *")
        sub.run(["conda activate phonopy"], shell=True)
        sub.run(["phonopy --qe -d --dim='{nd1} {nd2} {nd3}' -c {mat}.scf.in".format(nd1=prms.ndim[0], nd2=prms.ndim[1], nd3 = prms.ndim[2], mat=prms.mat)], shell=True)
        sub.run(["ls -lrt supercell-* > lsout"], shell=True)
        force_command: str = "phonopy -f "
        data:str = np.loadtxt("lsout",dtype="str",unpack=True,ndmin=0)
        files:str = data[len(data)-1]
        for fi in files:
            ncsc:str = fi.replace("supercell-","").replace(".in","")
            p = pathlib.Path(ncsc)
            if ( not p.is_dir() ):
                sub.run(["mkdir "+ncsc], shell=True)
                os.chdir(ncsc)
                sub.run(["cat ../header.in ../"+fi+" >| "+prms.mat+"-"+ncsc+".scf.in"], shell=True)
                # run scf calculation
                sub.run(["mpirun -np "+str(prms.nc)+" "+prms.exe+"/pw.x < "+prms.mat+"-"+ncsc+".scf.in > "+prms.mat+"-"+ncsc+".scf.out"], shell=True)
                sub.run(["rm work/"+prms.mat+".save/wfc*"], shell=True)
                force_command += ncsc+"/"+prms.mat+"-"+ncsc+".scf.out "
                os.chdir("../")
            
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
        self.dynmat: complex = []
        dynmat_data = data['phonon'][0]['dynamical_matrix']
        for row in dynmat_data:
            vals = np.reshape(row, (-1,2))
            self.dynmat.append(vals[:,0] + vals[:,1] * 1j)
        self.dynmat = np.array(self.dynmat)
        self.eigvals, self.eigvecs = np.linalg.eigh(self.dynmat)
        self.eigvecs = np.reshape(self.eigvecs, (3,len(self.eigvecs[0])//3,len(self.eigvecs)))
        self.freq: float = np.sqrt(np.abs(self.eigvals.real)) * np.sign(self.eigvals.real)
        self.freq = self.freq * const.conversion_factor_to_THz * const.THz2eV

        self.nmode: int = len(self.freq)
        self.qk: float = np.zeros(self.nmode)
        self.qk_force: float = np.zeros(self.nmode)
        (alat, plat, elements, nelems, natm, Rpos_g, volume) = subs.get_POSCAR("POSCAR_"+prms.stateg)
        (alat, plat, elements, nelems, natm, Rpos_e, volume) = subs.get_POSCAR("POSCAR_"+prms.statee)
        mass: float = []
        for i, nel in enumerate(nelems):
            for j in range(nel):
                mass.append(const.ELEMS_MASS[elements[i]]*const.uatm/const.me)
        if ( sw_qk in {"force","both"} ):
            Force_g: float = get_FORCE("FORCE_"+prms.stateg+".dat")
            Force_e: float = get_FORCE("FORCE_"+prms.statee+".dat")
        
        for i in range(natm):
            for j in range(3):
                if ( prms.sw_qk in {"pos","both"} ):
                    self.qk[:] += np.sqrt(mass[i]) * (Rpos_e[i,j]-Rpos_g[i,j]) * self.eigvecs[j,i,:].real / const.Bohr
                if ( prms.sw_qk in {"force","both"} ):
                    self.qk_force[:] += (1./((self.freq[:]/const.Ry)**2.)*np.sqrt(mass[i])) * (Force_e[i,j]-Force_g[i,j]) * self.eigvecs[j,i,:].real
        self.Sk: float = np.zeros(self.nmode)
        self.Sk_force: float = np.zeros(self.nmode)
        for i in range(self.nmode):
            self.Sk[i]: float = self.freq[i] * self.qk[i]**2. / ( 2.*const.Ry )
            self.Sk_force[i]: float = self.freq[i] * self.qk_force[i]**2. / ( 2.*const.Ry )
        energy: float = np.linspace(prms.emin_ph, prms.emax_ph, prms.ndive_ph)
        self.Sspec: float = np.zeros(prms.ndive)
        self.Sspec_force: float = np.zeros(prms.ndive)
        self.Stot: float = 0.0
        self.Stot_force: float = 0.0
        for i in range(self.nmode):
            if ( prms.sw_qk in {"pos","both"} ):
                self.Sspec[:] += self.Sk[i] * subs.Gaussian(energy[:]-self.freq[i])
                self.Stot += self.Sk[i]
            if ( prms.sw_qk in {"force","both"} ):
                self.Sspec_force[:] += self.Sk_force[i] * subs.Gaussian(energy[:]-self.freq[i])
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
