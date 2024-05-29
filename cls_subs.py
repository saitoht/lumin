import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import sys, os, pathlib

class subs:
    """ Class: define common subroutines """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor: initialize parameters"""

        """ Setting for plot """
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
        
        """ Physical constatns """
        self.Bohr = 0.529177210903        # \AA
        self.Ry = 13.6057039763           # eV
        self.uatm = 1.66054e-27           # kg
        self.me = 9.1093837e-31           # kg
        self.hconst = 4.135667696e-15     # eV s
        self.kB = 8.617333262e-5          # eV K^-1
        self.NA = 6.02214e23              # /mol
        self.hbar = 1.05457182e-34        # J s
        self.hbareV = 6.582119569e-16     # eV s
        self.ep = 1.60217663e-19          # C
        self.AU2GPa = 14710.513242194795  # A.U. to GPa
        self.eV2cm = 8065.73              # cm^-1
        self.THz2eV = 0.00413566553853599 # THz to eV
        self.fc2THz = 15.633302           # factor to convert into THz
        self.fceV2nm = 1239.8             # factor to convert eV into nm
        
        """ element mass """
        self.ELEMS_MASS = {"H":  1.00798, "He": 4.0026,  "Li": 6.968,   "Be": 9.01218, "B":  10.814,
                           "C":  12.0106, "N":  14.0069, "O":  15.9994, "F":  18.9984, "Ne": 20.1797,
                           "Na": 22.9898, "Mg": 24.306,  "Al": 26.9815, "Si": 28.085,  "P":  30.9738,
                           "S":  32.068,  "Cl": 35.452,  "Ar": 39.948,  "K":  39.0983, "Ca": 40.078,
                           "Sc": 44.9559, "Ti": 47.867,  "V": 50.9415,  "Cr": 51.9961, "Mn": 54.938,
                           "Fe": 55.845,  "Co": 58.9332, "Ni": 58.6934, "Cu": 63.546,  "Zn": 65.38,
                           "Ga": 69.723,  "Ge": 72.630,  "As": 74.9216, "Se": 78.971,  "Br": 79.904,
                           "Kr": 83.798,  "Rb": 85.4678, "Sr": 87.62,   "Y":  88.9058, "Zr": 91.224,
                           "Nb": 92.9064, "Mo": 95.95,   "Tc": 99.,     "Ru": 101.07,  "Rh": 102.906,
                           "Pd": 106.42,  "Ag": 107.868, "Cd": 112.414, "In": 114.818, "Sn": 118.710,
                           "Sb": 121.760, "Te": 127.60,  "I":  126.904, "Xe": 131.293, "Cs": 132.905,
                           "Ba": 137.327, "La": 138.905, "Ce": 140.116, "Pr": 140.908, "Nd": 144.242,
                           "Pm": 145,     "Sm": 150.36,  "Eu": 151.964, "Gd": 157.25,  "Tb": 158.925,
                           "Dy": 162.500, "Ho": 164.930, "Er": 167.259, "Tm": 168.934, "Yb": 173.045,
                           "Lu": 174.967, "Hf": 178.49,  "Ta": 180.948, "W":  183.84,  "Re": 186.207,
                           "Os": 190.23,  "Ir": 192.217, "Pt": 195.084, "Au": 196.967, "Hg": 200.592,
                           "Tl": 204.384, "Pb": 207.384, "Bi": 208.980, "Po": 210.,    "At": 210.,
                           "Rn": 222.,    "Fr": 223.,    "Ra": 226.,    "Ac": 227.,    "Th": 232.038,
                           "Pa": 231.038, "U" : 238.029, "Np": 237.,    "Pu": 239.,    "Am": 243.,
                           "Cm": 247.,    "Bk": 247.,    "Cf": 252.,    "Es": 252.,    "Fm": 257.,
                           "Md": 258.,    "No": 259.,    "Lr": 262.,    "Rf": 267.,    "Db": 268.,
                           "Sg": 271.,    "Bh": 272.,    "Hs": 277.,    "Mt": 276.,    "Ds": 281.,
                           "Rg": 280.,    "Cn": 285.,    "Nh": 278.,    "Fl": 289.,    "Mc": 289.,
                           "Lv": 293.,    "Ts": 293.,    "Og": 194.}

        """ Z number """
        self.Znum = {"H":  1,   "He": 2,   "Li": 3,   "Be": 4,   "B":  5,
                     "C":  6,   "N":  7,   "O":  8,   "F":  9,   "Ne": 10,
                     "Na": 11,  "Mg": 12,  "Al": 13,  "Si": 14,  "P":  15,
                     "S":  16,  "Cl": 17,  "Ar": 18,  "K":  19,  "Ca": 20,
                     "Sc": 21,  "Ti": 22,  "V": 23,   "Cr": 24,  "Mn": 25,
                     "Fe": 26,  "Co": 27,  "Ni": 28,  "Cu": 29,  "Zn": 30,
                     "Ga": 31,  "Ge": 32,  "As": 33,  "Se": 34,  "Br": 35,
                     "Kr": 36,  "Rb": 37,  "Sr": 38,  "Y": 39,   "Zr": 40,
                     "Nb": 41,  "Mo": 42,  "Tc": 43,  "Ru": 44,  "Rh": 45,
                     "Pd": 46,  "Ag": 47,  "Cd": 48,  "In": 49,  "Sn": 50,
                     "Sb": 51,  "Te": 52,  "I":  53,  "Xe": 54,  "Cs": 55,
                     "Ba": 56,  "La": 57,  "Ce": 58,  "Pr": 59,  "Nd": 60,
                     "Pm": 61,  "Sm": 62,  "Eu": 63,  "Gd": 64,  "Tb": 65,
                     "Dy": 66,  "Ho": 67,  "Er": 68,  "Tm": 69,  "Yb": 70,
                     "Lu": 71,  "Hf": 72,  "Ta": 73,  "W":  74,  "Re": 75,
                     "Os": 76,  "Ir": 77,  "Pt": 78,  "Au": 79,  "Hg": 80,
                     "Tl": 81,  "Pb": 82,  "Bi": 83,  "Po": 84,  "At": 85,
                     "Rn": 86,  "Fr": 87,  "Ra": 88,  "Ac": 89,  "Th": 90,
                     "Pa": 91,  "U" : 92,  "Np": 93,  "Pu": 94,  "Am": 95,
                     "Cm": 96,  "Bk": 97,  "Cf": 98,  "Es": 99,  "Fm": 100,
                     "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
                     "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
                     "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115,
                     "Lv": 116, "Ts": 117, "Og": 118}

        """ default values """
        """ for system """
        self.mat = ""
        self.nc = 20
        self.exe = "$HOME/qe/bin"
        self.code = "qe"
        self.nspin = 1
        """ for emod """
        self.brav = "cub"
        self.sw_plt_emod = False
        self.sw_Bmod = False
        self.sw_relax = False
        self.sw_egap = False
        self.sw_GW = False
        self.sw_auto = False
        self.dratio = 0.02
        self.ndiv_emod = 15
        """ for ccd """
        self.gamma = 1.0e-5
        self.Eabs0 = 2.146
        self.Eem0 = 1.702
        self.EFCg = 0.214
        self.dQ = 1.0
        self.include = None
        self.fac_include = 1.0
        self.sw_plt_ccd = False
        self.sw_eg = False
        self.sw_anharm = False
        self.sw_unit = "eV"
        self.sw_2dim = False
        self.sw_fix = "vertex"
        self.curvature = None
        self.dim_fit = 3
        self.stateg = "4A2g"
        self.statee = "4T2g"
        self.emin_ccd = 1.0
        self.emax_ccd = 3.0
        self.tempmin = 1.0
        self.tempmax = 1000.0
        self.temp = 298.0
        self.sw_spec_temp = False
        self.I0 = 1.0
        self.nmax = 30
        self.ndiv_e = 1000000
        self.ndiv_temp = 1000
        self.ndiv_eg = 12
        """ for ph """
        self.ndim = [1,1,1]
        self.sw_HR = "none"
        self.sw_phrun = False
        self.sw_plt_ph = False
        self.sigma = 6.0e-3
        self.ndiv_ph = 1000000
        self.emin_ph = 0.0
        self.emax_ph = 1.0
        self.tinf = 1.0e3
        self.ndiv_t = 1000000
        self.gamma_spec = 1.0e-5
        """ switch for execution """
        self.sw_run_emod = False
        self.sw_run_ccd = False
        self.sw_run_ph = False
        
        """ get parameters """
        subs.get_prms(self)
        
    ### ----------------------------------------------------------------------------- ###
    def get_prms(self):
        """ get parameters for calculations """
        
        sw_sys = False
        sw_emod = False
        sw_ph = False
        sw_ccd = False
        p = pathlib.Path("LUMIN")
        if ( p.exists() ):
            with open("LUMIN","r") as f:
                inputs = f.readlines()
            for ip in inputs:
                line = [ipp.replace("\n","") for ipp in ip.split()]
                if ( line[0] == "&sys" ):
                    sw_sys = True
                    sw_emod = False
                    sw_ph = False
                    sw_ccd = False
                elif ( line[0] == "&emod" ):
                    sw_sys = False
                    sw_emod = True
                    sw_ph = False
                    sw_ccd = False
                elif ( line[0] == "&ph" ):
                    sw_sys = False
                    sw_emod = False
                    sw_ph = True
                    sw_ccd = False
                elif ( line[0] == "&ccd" ):
                    sw_sys = False
                    sw_emod = False
                    sw_ph = False
                    sw_ccd = True
                else:
                    pass
                
                """ system """
                if ( sw_sys and line[0][0:4] == "mat=" ):
                    if ( len(line) == 2 ):
                        self.mat = line[1]
                    else:
                        self.mat = line[0].replace("mat=","")
                elif ( sw_sys and line[0][0:3] == "nc=" ):
                    if ( len(line) == 2 ):
                        self.nc = int(line[1])
                    else:
                        self.nc = int(line[0].replace("nc=",""))
                elif ( sw_sys and line[0][0:4] == "exe=" ):
                    if ( len(line) == 2 ):
                        self.exe = line[1]
                    else:
                        self.exe = line[0].replace("exe=","")
                elif ( sw_sys and line[0][0:5] == "code=" ):
                    if ( len(line) == 2 ):
                        self.code = line[1]
                    else:
                        self.code = line[0].replace("code=","")
                elif ( sw_sys and line[0][0:6] == "nspin=" ):
                    if ( len(line) == 2 ):
                        self.nspin = int(line[1])
                    else:
                        self.nspin = int(line[0].replace("nspin=",""))
                else:
                    pass

                """ elastic moduli """
                if ( sw_emod and line[0][0:5] == "brav=" ):
                    if ( len(line) == 2 ):
                        self.brav = line[1]
                    else:
                        self.brav = line[0].replace("brav=","")
                elif ( sw_emod and line[0][0:12] == "sw_plt_emod=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_emod = bool(int(line[1]))
                    else:
                        self.sw_plt_emod = bool(int(line[0].replace("sw_plt_emod=","")))
                elif ( sw_emod and line[0][0:8] == "sw_Bmod=" ):
                    if ( len(line) == 2 ):
                        self.sw_Bmod = bool(int(line[1]))
                    else:
                        self.sw_Bmod = bool(int(line[0].replace("sw_Bmod=","")))
                elif ( sw_emod and line[0][0:9] == "sw_relax=" ):
                    if ( len(line) == 2 ):
                        self.sw_relax = bool(int(line[1]))
                    else:
                        self.sw_relax = bool(int(line[0].replace("sw_relax=","")))
                elif ( sw_emod and line[0][0:5] == "sw_GW=" ):
                    if ( len(line) == 2 ):
                        self.sw_GW = bool(int(line[1]))
                    else:
                        self.sw_GW = bool(int(line[0].replace("sw_GW=","")))
                elif ( sw_emod and line[0][0:8] == "sw_auto=" ):
                    if ( len(line) == 2 ):
                        self.sw_auto = bool(int(line[1]))
                    else:
                        self.sw_auto = bool(int(line[0].replace("sw_auto=","")))
                elif ( sw_emod and line[0][0:8] == "sw_egap=" ):
                    if ( len(line) == 2 ):
                        self.sw_egap = bool(int(line[1]))
                    else:
                        self.sw_egap = bool(int(line[0].replace("sw_egap=","")))
                elif ( sw_emod and line[0][0:7] == "dratio=" ):
                    if ( len(line) == 2 ):
                        self.dratio = float(line[1])
                    else:
                        self.dratio = float(line[0].replace("dratio=",""))
                elif ( sw_emod and line[0][0:10] == "ndiv_emod=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_emod = int(line[1])
                    else:
                        self.ndiv_emod = int(line[0].replace("ndiv_emod=",""))
                else:
                    pass

                """ configuration coordinate diagram """
                if ( sw_ccd and line[0][0:6] == "Eabs0=" ):
                    if ( len(line) == 2 ):
                        self.Eabs0 = float(line[1])
                    else:
                        self.Eabs0 = float(line[0].replace("Eabs0=",""))
                elif ( sw_ccd and line[0][0:5] == "Eem0=" ):
                    if ( len(line) == 2 ):
                        self.Eem0 = float(line[1])
                    else:
                        self.Eem0 = float(line[0].replace("Eem0=",""))
                elif ( sw_ccd and line[0][0:5] == "EFCg=" ):
                    if ( len(line) == 2 ):
                        self.EFCg = float(line[1])
                    else:
                        self.EFCg = float(line[0].replace("EFCg=",""))
                elif ( sw_ccd and line[0][0:3] == "dQ=" ):
                    if ( len(line) == 2 ):
                        self.dQ = float(line[1])
                    else:
                        self.dQ = float(line[0].replace("dQ=",""))
                elif ( sw_ccd and line[0][0:3] == "I0=" ):
                    if ( len(line) == 2 ):
                        self.I0 = float(line[1])
                    else:
                        self.I0 = float(line[0].replace("I0=",""))
                elif ( sw_ccd and line[0][0:6] == "gamma=" ):
                    if ( len(line) == 2 ):
                        self.gamma = float(line[1])
                    else:
                        self.gamma = float(line[0].replace("gamma=",""))
                elif ( sw_ccd and line[0][0:11] == "sw_plt_ccd=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_ccd = bool(int(line[1]))
                    else:
                        self.sw_plt_ccd = bool(int(line[0].replace("sw_plt_ccd=","")))
                elif ( sw_ccd and line[0][0:6] == "sw_eg=" ):
                    if ( len(line) == 2 ):
                        self.sw_eg = bool(int(line[1]))
                    else:
                        self.sw_eg = bool(int(line[0].replace("sw_eg=","")))
                elif ( sw_ccd and line[0][0:10] == "sw_anharm=" ):
                    if ( len(line) == 2 ):
                        self.sw_anharm = bool(int(line[1]))
                    else:
                        self.sw_anharm = bool(int(line[0].replace("sw_anharm=","")))
                elif ( sw_ccd and line[0][0:13] == "sw_spec_temp=" ):
                    if ( len(line) == 2 ):
                        self.sw_spec_temp = bool(int(line[1]))
                    else:
                        self.sw_spec_temp = bool(int(line[0].replace("sw_spec_temp=","")))
                elif ( sw_ccd and line[0][0:8] == "sw_unit=" ):
                    if ( len(line) == 2 ):
                        self.sw_unit = line[1]
                    else:
                        self.sw_unit = line[0].replace("sw_unit=","")
                elif ( sw_ccd and line[0][0:8] == "include=" ):
                    if ( len(line) == 2 ):
                        self.include = line[1]
                    else:
                        self.include = line[0].replace("include=","")
                elif ( sw_ccd and line[0][0:12] == "fac_include=" ):
                    if ( len(line) == 2 ):
                        self.fac_include = float(line[1])
                    else:
                        self.fac_include = float(line[0].replace("fac_include=",""))
                elif ( sw_ccd and line[0][0:8] == "sw_2dim=" ):
                    if ( len(line) == 2 ):
                        self.sw_2dim = bool(int(line[1]))
                    else:
                        self.sw_2dim = bool(int(line[0].replace("sw_2dim=","")))
                elif ( sw_ccd and line[0][0:7] == "sw_fix=" ):
                    if ( len(line) == 2 ):
                        self.sw_fix = line[1]
                    else:
                        self.sw_fix = line[0].replace("sw_fix=","")
                elif ( sw_ccd and line[0][0:10] == "curvature=" ):
                    if ( len(line) == 2 ):
                        self.curvature = float(line[1])
                    else:
                        self.curvature = float(line[0].replace("curvature=",""))
                elif ( sw_ccd and line[0][0:8] == "dim_fit=" ):
                    if ( len(line) == 2 ):
                        self.dim_fit = int(line[1])
                    else:
                        self.dim_fit = int(line[0].replace("dim_fit=",""))
                elif ( sw_ccd and line[0][0:7] == "stateg=" ):
                    if ( len(line) == 2 ):
                        self.stateg = line[1]
                    else:
                        self.stateg = line[0].replace("stateg=","")
                elif ( sw_ccd and line[0][0:7] == "statee=" ):
                    if ( len(line) == 2 ):
                        self.statee = line[1]
                    else:
                        self.statee = line[0].replace("statee=","")
                elif ( sw_ccd and line[0][0:9] == "emin_ccd=" ):
                    if ( len(line) == 2 ):
                        self.emin_ccd = float(line[1])
                    else:
                        self.emin_ccd = float(line[0].replace("emin_ccd=",""))
                elif ( sw_ccd and line[0][0:9] == "emax_ccd=" ):
                    if ( len(line) == 2 ):
                        self.emax_ccd = float(line[1])
                    else:
                        self.emax_ccd = float(line[0].replace("emax_ccd=",""))
                elif ( sw_ccd and line[0][0:8] == "tempmin=" ):
                    if ( len(line) == 2 ):
                        self.tempmin = float(line[1])
                    else:
                        self.tempmin = float(line[0].replace("tempmin=",""))
                elif ( sw_ccd and line[0][0:8] == "tempmax=" ):
                    if ( len(line) == 2 ):
                        self.tempmax = float(line[1])
                    else:
                        self.tempmax = float(line[0].replace("tempmax=",""))
                elif ( sw_ccd and line[0][0:8] == "temp=" ):
                    if ( len(line) == 2 ):
                        self.temp = float(line[1])
                    else:
                        self.temp = float(line[0].replace("temp=",""))
                elif ( sw_ccd and line[0][0:5] == "nmax=" ):
                    if ( len(line) == 2 ):
                        self.nmax = int(line[1])
                    else:
                        self.nmax = int(line[0].replace("nmax=",""))
                elif ( sw_ccd and line[0][0:7] == "ndiv_e=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_e = int(line[1])
                    else:
                        self.ndiv_e = int(line[0].replace("ndiv_e=",""))
                elif ( sw_ccd and line[0][0:10] == "ndiv_temp=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_temp = int(line[1])
                    else:
                        self.ndiv_temp = int(line[0].replace("ndiv_temp=",""))
                elif ( sw_ccd and line[0][0:8] == "ndiv_eg=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_eg = int(line[1])
                    else:
                        self.ndiv_eg = int(line[0].replace("ndiv_eg=",""))
                else:
                    pass
            
                """ phonon """
                if ( sw_ph and line[0][0:5] == "ndim=" ):
                    if ( len(line) == 4 ):
                        self.ndim = [int(l) for l in line[1:4]]
                    else:
                        self.ndim = [int(line[0].replace("ndim=","")),int(line[1]),int(line[2])]
                    self.ndim = np.array(self.ndim)
                elif ( sw_ph and line[0][0:7] == "emin_ph=" ):
                    if ( len(line) == 2 ):
                        self.emin_ph = float(line[1])
                    else:
                        self.emin_ph = float(line[0].replace("emin_ph=",""))
                elif ( sw_ph and line[0][0:7] == "emax_ph=" ):
                    if ( len(line) == 2 ):
                        self.emax_ph = float(line[1])
                    else:
                        self.emax_ph = float(line[0].replace("emax_ph=",""))
                elif ( sw_ph and line[0][0:6] == "sigma=" ):
                    if ( len(line) == 2 ):
                        self.sigma = float(line[1])
                    else:
                        self.sigma = float(line[0].replace("sigma=",""))
                elif ( sw_ph and line[0][0:5] == "tinf=" ):
                    if ( len(line) == 2 ):
                        self.tinf = float(line[1])
                    else:
                        self.tinf = float(line[0].replace("tinf=",""))
                elif ( sw_ph and line[0][0:6] == "ndiv_t=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_t = int(line[1])
                    else:
                        self.ndiv_t = int(line[0].replace("ndiv_t=",""))
                elif ( sw_ph and line[0][0:11] == "gamma_spec=" ):
                    if ( len(line) == 2 ):
                        self.gamma_spec = float(line[1])
                    else:
                        self.gamma_spec = float(line[0].replace("gamma_spec=",""))
                elif ( sw_ph and line[0][0:6] == "sw_HR=" ):
                    if ( len(line) == 2 ):
                        self.sw_HR = line[1]
                    else:
                        self.sw_HR = line[0].replace("sw_HR=","")
                elif ( sw_ph and line[0][0:9] == "sw_phrun=" ):
                    if ( len(line) == 2 ):
                        self.sw_phrun = bool(int(line[1]))
                    else:
                        self.sw_phrun = bool(int(line[0].replace("sw_phrun=","")))
                elif ( sw_ph and line[0][0:10] == "sw_plt_ph=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_ph = line[1]
                    else:
                        self.sw_plt_ph = line[0].replace("sw_plt_ph=","")
                else:
                    pass
        else:
            pass

        args = sys.argv
        if ( (len(args) == 1 and self.mat == "") or "-h" in args ):
            helpmessage: str = """ *** USAGE: python lumin.py {mat} [-h] [-nc {nc}] [-exe {exe}] [-p] [-emod] [-ccd] [-ph]
                                -h : print helpmessage
                                -nc {nc} : number of cpu {20}
                                -exe {exe} : execute path {$HOME/qe/bin}
                                -code {code} : code for DFT (qe, ecalj) {qe}
                                -nsp {nspin} : number of spin states under consideration {1}
                                -p : sw_plt == True for all script {False}
                        [-emod] -pemod : sw_plt_emod == True {False}
                                -brav {brav} : Bravias type (cub, tet, hex, mono, tri, ...) {cub}
                                -dr {dratio} : distortion ratio {0.02}
                                -ndem {ndiv_emod} : number of devision for elastic moduli calc. {15}
                                -Bmod : calculate Bmod or not {False}
                                -relax : consider relaxation of atomic positions or not {False}
                                -egap : extract band gap energy {False}
                                -GW : perform GW calculation, only valid for code=ecalj case {False}
                                -auto : automatic evaluation of elastic moduli & energy band gap {False}
                        [-ccd]  -pccd : sw_plt_ccd == True {False}
                                -Eabs {Eabs0} : absorption energy by first-principles calc. {2.146}
                                -Eem {Eem0} : emission energy by first-principles calc. {1.702}
                                -FC {EFCg} : Frank-Condon params. for ground state energy curve {0.214}
                                -dQ {dQ} : difference of the configuration coordinate b/w exited & ground for intemediate state calculations {1.0}
                                -I0 {I0} : Intensity of the spectrum {1.0}
                                -gamma {gamma} : broadness params. of Lorentzian {1.0e-5}
                                -eg : generate intermediate structures b/w excited & ground or not {False}
                                -anharm : consider anharmonicity in adiabatic energy curve or not {False}
                                -2dim : plot 2dim fit or not {False}
                                -fix {sw_fix} : switch for 2dim fit in 1D-CCD, fix vertex position (vertex) or curvature (curve) {vertex}
                                -curv {curvature} : curvature to be fixed {None}
                                -dim {dim_fit} : dimension to fit energy curve {3}
                                -stateg {stateg} : name of the ground state {4A2g}
                                -statee {statee} : name of the excited state {4T2g}
                                -include {include} : name of file to be shown in spectrum {None}
                                -fac {fac_include} : factor of include file x axis {1.0}
                                -unit {sw_unit} : select unit of quantities (eV, cm^-1, nm) {eV}
                                -emin_ccd {emin_ccd} : energy min. for plotting {1.0}
                                -emax_ccd {emax_ccd} : energy max. for plotting {3.0}
                                -tempmin {tempmin} : temperature min. {1.0}
                                -tempmax {tempmax} : temperature max. {1.0e3}
                                -temp {temp} : temperature {298.0}
                                -spec_temp : show spectrum at finite temp or not {False}
                                -nmax {nmax} : cut-off number of summation, larger value is plausible {30}
                                -nde {ndiv_e} : number of devision for energy {1000000}
                                -ndtmp {ndiv_temp} : number of devision for temperature {1000}
                                -ndeg {ndiv_eg} : number of devision for intermediate structures {12}
                        [-ph]   -pph : sw_plt_ph == True {False}
                                -ndim {ndim[1:3]} : number of dimension for generating supercell {1 1 1}
                                -HR {sw_HR} : switch for means to evaluate q_k (none, pos, force, both) {none}
                                -phrun : perform phonon calc. or not {False}
                                -nde_ph {ndiv_ph} : number of division for phonon spectrum {1000000}
                                -emin_ph {emin_ph} : energy min. for phonon spectrum {0.0}
                                -emax_ph {emax_ph} : energy max. for phonon spectrum {1.0}
                                -sigma {sigma} : standard derivation of Gaussian {6.0e-3}
                                -tinf {tinf} : inf value of t {1.0e3}
                                -ndt {ndiv_t} : number of division for t {1000000}
                                -gsp {gamma_spec} : broadening parameter of spectrum {1.0e-5} """
            print(helpmessage)
            sys.exit()

        if ( ( self.mat == "" ) and ( not args[1][0] == "-" ) ):
            self.mat: str = args[1]
        if ( "-emod" in args ):
            self.sw_run_emod = True
        if ( "-ccd" in args ):
            self.sw_run_ccd = True
        if ( "-ph" in args ):
            self.sw_run_ph = True
        if ( "-nc" in args ):
            self.nc = int(args[args.index("-nc")+1])
        if ( "-exe" in args ):
            self.exe = args[args.index("-exe")+1]
        if ( "-code" in args ):
            self.code = args[args.index("-code")+1]
        if ( "-nsp" in args ):
            self.nspin = int(args[args.index("-nsp")+1])
        if ( "-p" in args ):
            self.sw_plt_emod = True
            self.sw_plt_ccd = True
            self.sw_plt_ph = True
        if ( "-pemod" in args ):
            self.sw_plt_emod = True
        if ( "-pccd" in args ):
            self.sw_plt_ccd = True
        if ( "-pph" in args ):
            self.sw_plt_ph = True
        """ emod """
        if ( "-brav" in args ):
            self.brav = args[args.index("-brav")+1]
        if ( "-dr" in args ):
            self.dratio = float(args[args.index("-dr")+1])
        if ( "-ndem" in args ):
            self.ndiv_emod = int(args[args.index("-ndem")+1])
        if ( "-Bmod" in args ):
            self.sw_Bmod = True
        if ( "-relax" in args ):
            self.sw_relax = True
        if ( "-egap" in args ):
            self.sw_egap = True
        if ( "-GW" in args ):
            self.sw_GW = True
        if ( "-auto" in args ):
            self.sw_auto = True
        """ ccd """
        if ( "-Eabs" in args ):
            self.Eabs0 = float(args[args.index("-Eabs")+1])
        if ( "-Eem" in args ):
            self.Eem0 = float(args[args.index("-Eem")+1])
        if ( "-FC" in args ):
            self.EFCg = float(args[args.index("-FC")+1])
        if ( "-dQ" in args ):
            self.dQ = float(args[args.index("-dQ")+1])
        if ( "-I0" in args ):
            self.I0 = float(args[args.index("-I0")+1])
        if ( "-gamma" in args ):
            self.gamma = float(args[args.index("-gamma")+1])
        if ( "-eg" in args ):
            self.sw_eg = True
        if ( "-anharm" in args ):
            self.sw_anharm = True
        if ( "-2dim" in args ):
            self.sw_2dim = True
        if ( "-spec_temp" in args ):
            self.sw_spec_temp = True
        if ( "-stateg" in args ):
            self.stateg = args[args.index("-stateg")+1]
        if ( "-statee" in args ):
            self.stateg = args[args.index("-statee")+1]
        if ( "-unit" in args ):
            self.sw_unit = args[args.index("-unit")+1]
        if ( "-fix" in args ):
            self.sw_fix = args[args.index("-fix")+1]
        if ( "-curv" in args ):
            self.curvature = float(args[args.index("-curv")+1])
        if ( "-include" in args ):
            self.include = args[args.index("-include")+1]
        if ( "-fac" in args ):
            self.fac_include = float(args[args.index("-fac")+1])
        if ( "-dim" in args ):
            self.dim_fit = int(args[args.index("-dim")+1])
        if ( "-emin_ccd" in args ):
            self.emin_ccd = float(args[args.index("-emin_ccd")+1])
        if ( "-emax_ccd" in args ):
            self.emax_ccd = float(args[args.index("-emax_ccd")+1])
        if ( "-tempmin" in args ):
            self.tempmin = float(args[args.index("-tempmin")+1])
        if ( "-tempmax" in args ):
            self.tempmax = float(args[args.index("-tempmax")+1])
        if ( "-temp" in args ):
            self.temp = float(args[args.index("-temp")+1])
        if ( "-nmax" in args ):
            self.nmax = int(args[args.index("-nmax")+1])
        if ( "-nde" in args ):
            self.ndiv_e = int(args[args.index("-nde")+1])
        if ( "-ndtmp" in args ):
            self.ndiv_temp = int(args[args.index("-ndtmp")+1])
        if ( "-ndeg" in args ):
            self.ndiv_eg = int(args[args.index("-ndeg")+1])
        """ ph """
        if ( "-ndim" in args ):
            self.ndim = [int(args[args.index("-ndim")+1+i]) for i in range(3)]
            self.ndim = np.array(self.ndim)
        if ( "-HR" in args ):
            self.sw_HR = args[args.index("-HR")+1]
        if ( "-phrun" in args ):
            self.sw_phrun = True
        if ( "-nde_ph" in args ):
            self.ndiv_ph = int(args[args.index("-nde_ph")+1])
        if ( "-emin_ph" in args ):
            self.emin_ph = float(args[args.index("-emin_ph")+1])
        if ( "-emax_ph" in args ):
            self.emax_ph = float(args[args.index("-emax_ph")+1])
        if ( "-sigma" in args ):
            self.sigma = float(args[args.index("-sigma")+1])
        if ( "-tinf" in args ):
            self.tinf = float(args[args.index("-tinf")+1])
        if ( "-gsp" in args ):
            self.gamma_spec = float(args[args.index("-gsp")+1])
        if ( "-ndt" in args ):
            self.ndiv_t = int(args[args.index("-ndt")+1])

    ### ----------------------------------------------------------------------------- ###
    def ck_prms(self):
        """ check parameters """

        print("* --- check parameters --- *")
        print("*** system parameters ***")
        print("* mat: ", self.mat)
        print("* sw_run_emod: ", self.sw_run_emod)
        print("* sw_run_ccd: ", self.sw_run_ccd)
        print("* sw_run_ph: ", self.sw_run_ph)
        print("* nc: ", self.nc)
        print("* exe: ", self.exe)
        print("* code: ", self.code)
        print("* nspin: ", self.nspin)
        print("*")
        if ( self.sw_run_emod ):
            print("*** elastic moduli parameters ***")
            print("* brav: ", self.brav)
            print("* sw_plt_emod: ", self.sw_plt_emod)
            print("* sw_Bmod: ", self.sw_Bmod)
            print("* sw_relax: ", self.sw_relax)
            print("* sw_egap: ", self.sw_egap)
            print("* sw_GW: ", self.sw_GW)
            print("* sw_auto: ", self.sw_auto)
            print("* dratio: ", self.dratio)
            print("* ndiv_emod: ", self.ndiv_emod)
            print("*")
        if ( self.sw_run_ccd ):
            print("*** configuration coordinate diagram parameters ***")        
            print("* Eabs0 (eV): ", self.Eabs0)
            print("* Eem0 (eV): ", self.Eem0)
            print("* EFCg (eV): ", self.EFCg)
            print("* dQ (eV): ", self.dQ)
            print("* I0: ", self.I0)
            print("* gamma: ", self.gamma)
            print("* sw_plt_ccd: ", self.sw_plt_ccd)
            print("* sw_eg: ", self.sw_eg)
            print("* sw_anharm: ", self.sw_anharm)
            print("* sw_2dim: ", self.sw_2dim)
            print("* sw_fix: ", self.sw_fix)
            print("* sw_unit: ", self.sw_unit)
            print("* curvature: ", self.curvature)
            print("* include: ", self.include)
            print("* dim_fit: ", self.dim_fit)
            print("* fac_include: ", self.fac_include)
            print("* stateg: ", self.stateg)
            print("* statee: ", self.statee)
            print("* emin_ccd (eV): ", self.emin_ccd)
            print("* emax_ccd (eV): ", self.emax_ccd)
            print("* tempmin (K): ", self.tempmin)
            print("* tempmax (K): ", self.tempmax)
            print("* temp (K): ", self.temp)
            print("* sw_spec_temp: ", self.sw_spec_temp)
            print("* nmax: ", self.nmax)
            print("* ndiv_e: ", self.ndiv_e)
            print("* ndiv_temp: ", self.ndiv_temp)
            print("* ndiv_eg: ", self.ndiv_eg)
            print("*")
        if ( self.sw_run_ph ):
            print("*** phonon parameters ***")
            print("* ndim: ", self.ndim)
            print("* emin_ph (eV): ", self.emin_ph)
            print("* emax_ph (eV): ", self.emax_ph)
            print("* sigma: ", self.sigma)
            print("* tinf: ", self.tinf)
            print("* ndiv_t: ", self.ndiv_t)
            print("* gamma_spec: ", self.gamma_spec)
            print("* sw_HR: ", self.sw_HR)
            print("* sw_phrun: ", self.sw_phrun)
            print("*")
        print("* --- Finish check parameters --- *")

    ### ----------------------------------------------------------------------------- ###
    def set_prms(self, material, bravias):
        """ set parameters """
    
        self.mat = material
        self.brav = bravias

    ### ----------------------------------------------------------------------------- ###
    def set_prms_ccd(self, Eabs, Eem, EFC):
        """ set ccd parameters """

        self.Eabs0 = Eabs
        self.Eem0 = Eem
        self.EFCg = EFC
        
    ### ----------------------------------------------------------------------------- ###
    def get_POSCAR(self, fn):
        """ get information from POSCAR """
        
        with open(fn, "r") as f:
            data = f.readlines()
        alat = float(data[1].split()[0]) / self.Bohr
        plat = [[float(d) for d in data[2].split()],
                      [float(d) for d in data[3].split()],
                      [float(d) for d in data[4].split()]]
        elements = [d for d in data[5].split()]
        nelems = [int(d) for d in data[6].split()]
        natm = sum(nelems)
        pos = []
        for dat in data[8:]:
            pos.append([float(d) for d in dat.split()])
        plat = np.array(plat)
        pos = np.array(pos)
        volume = (alat**3.) * abs(np.dot(plat[0],np.cross(plat[1],plat[2])))
        return (alat, plat, elements, nelems, natm, pos, volume)

    ### ----------------------------------------------------------------------------- ###
    def get_FORCE(self, fn):
        """ get force information of ground & excited states """
        
        p = pathlib.Path(fn)
        if ( not p.exists() ):
            print("*** ERROR: {fn} does not exist!!!".format(fn=fn))
            sys.exit()
        data: str = np.loadtxt(fn,dtype="str",unpack=True,ndmin=0)
        Force: float = [[float(data[6,i]),float(data[7,i]),float(data[8,i])] for i in range(len(data[0]))]
        Force = np.array(Force)
        return Force

    ### ----------------------------------------------------------------------------- ###
    def Lorentzian(self, x):
        """ define Lorentzian """
        
        return (1./np.pi)*(0.5*self.gamma / (x**2.+(0.5*self.gamma)**2.))

    ### ----------------------------------------------------------------------------- ###
    def Gaussian(self, x):
        """ define Gaussian """
        
        return (1./(self.sigma*np.sqrt(2.*np.pi))) * np.exp(-x**2./(2.*self.sigma**2.))

    ### ----------------------------------------------------------------------------- ###
    def E2lambda(self, energy):
        """ transfrom energy to wavelength """
        
        return self.fceV2nm/energy  # nm

    ### ----------------------------------------------------------------------------- ###
    def lambda2E(self, lam):
        """ transform wavelength to energy """
        
        return self.fceV2nm/lam     # eV
