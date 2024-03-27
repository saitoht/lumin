import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import sys, os, pathlib

class subs:
    """ Class: define common subroutines """
    ### ----------------------------------------------------------------------------- ###
    def __init__(self):
        """ Constructor: initialize parameters"""
        
        """ Physical constatns """
        self.Bohr:float = 0.529177210903        # \AA
        self.Ry:float = 13.6057039763           # eV
        self.uatm:float = 1.66054e-27           # kg
        self.me:float = 9.1093837e-31           # kg
        self.hconst:float = 4.135667696e-15     # eV s
        self.kB:float = 8.617333262e-5          # eV K^-1
        self.NA:float = 6.02214e23              # /mol
        self.hbar:float = 1.05457182e-34        # J s
        self.ep:float = 1.60217663e-19          # C
        self.AU2GPa:float = 14710.513242194795  # A.U. to GPa
        self.eV2cm:float = 8065.73              # cm^-1
        self.THz2eV:float = 0.00413566553853599 # THz to eV
        self.fc2THz:float = 15.633302           # factor to convert into THz
        self.fceV2nm:float = 1239.8             # factor to convert eV into nm
        
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

        """ default values """
        """ for system """
        self.mat:int = ""
        self.nc:int = 20
        self.exe:int = "$HOME/qe/bin"
        """ for emod """
        self.brav:str = "cub"
        self.sw_plt_emod:bool = False
        self.sw_Bmod:bool = False
        self.sw_egap:bool = False
        self.dratio:float = 0.02
        self.ndiv_emod:int = 15
        """ for ccd """
        self.gamma:float = 1.0e-5
        self.Eabs0:float = 2.146
        self.Eem0:float = 1.702
        self.EFCg:float = 0.214
        self.dQ:float = 1.0
        self.sw_plt_ccd:bool = False
        self.sw_eg:bool = False
        self.sw_unit:str = "eV"
        self.stateg:str = "4A2g"
        self.statee:str = "4T2g"
        self.emin_ccd:float = 1.0
        self.emax_ccd:float = 3.0
        self.tempmin:float = 1.0
        self.tempmax:float = 1000
        self.I0:float = 1.0
        self.nmax:int = 30
        self.ndiv_e:int = 1000000
        self.ndiv_temp:int = 1000
        self.ndiv_eg:int = 12
        """ for ph """
        self.ndim:int = [1,1,1]
        self.sw_HR:str = "none"
        self.sw_phrun:bool = False
        self.sw_plt_ph:bool = False
        self.sigma:float = 6.0e-3
        self.ndiv_ph:int = 1000000
        self.emin_ph:float = 0.0
        self.emax_ph:float = 1.0
        """ switch for execution """
        self.sw_run_emod:bool = False
        self.sw_run_ccd:bool = False
        self.sw_run_ph:bool = False
        
        """ get parameters """
        subs.get_prms(self)
        
    ### ----------------------------------------------------------------------------- ###
    def get_prms(self):
        """ get parameters for calculations """
        
        sw_sys:bool = False
        sw_emod:bool = False
        sw_ph:bool = False
        sw_ccd:bool = False
        p = pathlib.Path("LUMIN")
        if ( p.exists() ):
            with open("LUMIN","r") as f:
                inputs:str = f.readlines()
            for ip in inputs:
                line:str = [ipp.replace("\n","") for ipp in ip.split()]
                if ( line[0] == "&sys" ):
                    sw_sys:bool = True
                    sw_emod:bool = False
                    sw_ph:bool = False
                    sw_ccd:bool = False
                elif ( line[0] == "&emod" ):
                    sw_sys:bool = False
                    sw_emod:bool = True
                    sw_ph:bool = False
                    sw_ccd:bool = False
                elif ( line[0] == "&ph" ):
                    sw_sys:bool = False
                    sw_emod:bool = False
                    sw_ph:bool = True
                    sw_ccd:bool = False
                elif ( line[0] == "&ccd" ):
                    sw_sys:bool = False
                    sw_emod:bool = False
                    sw_ph:bool = False
                    sw_ccd:bool = True
                else:
                    pass
                
                """ system """
                if ( sw_sys and line[0][0:4] == "mat=" ):
                    if ( len(line) == 2 ):
                        self.mat:str = line[1]
                    else:
                        self.mat:str = line[0].replace("mat=","")
                elif ( sw_sys and line[0][0:3] == "nc=" ):
                    if ( len(line) == 2 ):
                        self.nc:int = int(line[1])
                    else:
                        self.nc:int = int(line[0].replace("nc=",""))
                elif ( sw_sys and line[0][0:4] == "exe=" ):
                    if ( len(line) == 2 ):
                        self.exe:str = line[1]
                    else:
                        self.exe:str = line[0].replace("exe=","")
                else:
                    pass

                """ elastic moduli """
                if ( sw_emod and line[0][0:5] == "brav=" ):
                    if ( len(line) == 2 ):
                        self.brav:str = line[1]
                    else:
                        self.brav:str = line[0].replace("brav=","")
                elif ( sw_emod and line[0][0:12] == "sw_plt_emod=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_emod:bool = bool(int(line[1]))
                    else:
                        self.sw_plt_emod:bool = bool(int(line[0].replace("sw_plt_emod=","")))
                elif ( sw_emod and line[0][0:8] == "sw_Bmod=" ):
                    if ( len(line) == 2 ):
                        self.sw_Bmod:bool = bool(int(line[1]))
                    else:
                        self.sw_Bmod:bool = bool(int(line[0].replace("sw_Bmod=","")))
                elif ( sw_emod and line[0][0:8] == "sw_egap=" ):
                    if ( len(line) == 2 ):
                        self.sw_egap:bool = bool(int(line[1]))
                    else:
                        self.sw_egap:bool = bool(int(line[0].replace("sw_egap=","")))
                elif ( sw_emod and line[0][0:7] == "dratio=" ):
                    if ( len(line) == 2 ):
                        self.dratio:float = float(line[1])
                    else:
                        self.dratio:float = float(line[0].replace("dratio=",""))
                elif ( sw_emod and line[0][0:10] == "ndiv_emod=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_emod:int = int(line[1])
                    else:
                        self.ndiv_emod:int = int(line[0].replace("ndiv_emod=",""))
                else:
                    pass

                """ configuration coordinate diagram """
                if ( sw_ccd and line[0][0:6] == "Eabs0=" ):
                    if ( len(line) == 2 ):
                        self.Eabs0:float = float(line[1])
                    else:
                        self.Eabs0:float = float(line[0].replace("Eabs0=",""))
                elif ( sw_ccd and line[0][0:5] == "Eem0=" ):
                    if ( len(line) == 2 ):
                        self.Eem0:float = float(line[1])
                    else:
                        self.Eem0:float = float(line[0].replace("Eem0=",""))
                elif ( sw_ccd and line[0][0:5] == "EFCg=" ):
                    if ( len(line) == 2 ):
                        self.EFCg:float = float(line[1])
                    else:
                        self.EFCg:float = float(line[0].replace("EFCg=",""))
                elif ( sw_ccd and line[0][0:3] == "dQ=" ):
                    if ( len(line) == 2 ):
                        self.dQ:float = float(line[1])
                    else:
                        self.dQ:float = float(line[0].replace("dQ=",""))
                elif ( sw_ccd and line[0][0:3] == "I0=" ):
                    if ( len(line) == 2 ):
                        self.I0:float = float(line[1])
                    else:
                        self.I0:float = float(line[0].replace("I0=",""))
                elif ( sw_ccd and line[0][0:6] == "gamma=" ):
                    if ( len(line) == 2 ):
                        self.gamma:float = float(line[1])
                    else:
                        self.gamma:float = float(line[0].replace("gamma=",""))
                elif ( sw_ccd and line[0][0:11] == "sw_plt_ccd=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_ccd:bool = bool(int(line[1]))
                    else:
                        self.sw_plt_ccd:bool = bool(int(line[0].replace("sw_plt_ccd=","")))
                elif ( sw_ccd and line[0][0:6] == "sw_eg=" ):
                    if ( len(line) == 2 ):
                        self.sw_eg:bool = bool(int(line[1]))
                    else:
                        self.sw_eg:bool = bool(int(line[0].replace("sw_eg=","")))
                elif ( sw_ccd and line[0][0:8] == "sw_unit=" ):
                    if ( len(line) == 2 ):
                        self.sw_unit:str = line[1]
                    else:
                        self.sw_unit:str = line[0].replace("sw_unit=","")
                elif ( sw_ccd and line[0][0:7] == "stateg=" ):
                    if ( len(line) == 2 ):
                        self.stateg:str = line[1]
                    else:
                        self.stateg:str = line[0].replace("stateg=","")
                elif ( sw_ccd and line[0][0:7] == "statee=" ):
                    if ( len(line) == 2 ):
                        self.statee:str = line[1]
                    else:
                        self.statee:str = line[0].replace("statee=","")
                elif ( sw_ccd and line[0][0:9] == "emin_ccd=" ):
                    if ( len(line) == 2 ):
                        self.emin_ccd:float = float(line[1])
                    else:
                        self.emin_ccd:float = float(line[0].replace("emin_ccd=",""))
                elif ( sw_ccd and line[0][0:9] == "emax_ccd=" ):
                    if ( len(line) == 2 ):
                        self.emax_ccd:float = float(line[1])
                    else:
                        self.emax_ccd:float = float(line[0].replace("emax_ccd=",""))
                elif ( sw_ccd and line[0][0:8] == "tempmin=" ):
                    if ( len(line) == 2 ):
                        self.tempmin:float = float(line[1])
                    else:
                        self.tempmin:float = float(line[0].replace("tempmin=",""))
                elif ( sw_ccd and line[0][0:8] == "tempmax=" ):
                    if ( len(line) == 2 ):
                        self.tempmax:float = float(line[1])
                    else:
                        self.tempmax:float = float(line[0].replace("tempmax=",""))
                elif ( sw_ccd and line[0][0:5] == "nmax=" ):
                    if ( len(line) == 2 ):
                        self.nmax:int = int(line[1])
                    else:
                        self.nmax:int = int(line[0].replace("nmax=",""))
                elif ( sw_ccd and line[0][0:7] == "ndiv_e=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_e:int = int(line[1])
                    else:
                        self.ndiv_e:int = int(line[0].replace("ndiv_e=",""))
                elif ( sw_ccd and line[0][0:10] == "ndiv_temp=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_temp:int = int(line[1])
                    else:
                        self.ndiv_temp:int = int(line[0].replace("ndiv_temp=",""))
                elif ( sw_ccd and line[0][0:8] == "ndiv_eg=" ):
                    if ( len(line) == 2 ):
                        self.ndiv_eg:int = int(line[1])
                    else:
                        self.ndiv_eg:int = int(line[0].replace("ndiv_eg=",""))
                else:
                    pass
            
                """ phonon """
                if ( sw_ph and line[0][0:5] == "ndim=" ):
                    if ( len(line) == 4 ):
                        self.ndim:int = [int(l) for l in line[1:4]]
                    else:
                        self.ndim:int = [int(line[0].replace("ndim=","")),int(line[1]),int(line[2])]
                    self.ndim = np.array(self.ndim)
                elif ( sw_ph and line[0][0:7] == "emin_ph=" ):
                    if ( len(line) == 2 ):
                        self.emin_ph:float = float(line[1])
                    else:
                        self.emin_ph:float = float(line[0].replace("emin_ph=",""))
                elif ( sw_ph and line[0][0:7] == "emax_ph=" ):
                    if ( len(line) == 2 ):
                        self.emax_ph:float = float(line[1])
                    else:
                        self.emax_ph:float = float(line[0].replace("emax_ph=",""))
                elif ( sw_ph and line[0][0:6] == "sigma=" ):
                    if ( len(line) == 2 ):
                        self.sigma:float = float(line[1])
                    else:
                        self.sigma:float = float(line[0].replace("sigma=",""))
                elif ( sw_ph and line[0][0:6] == "sw_HR=" ):
                    if ( len(line) == 2 ):
                        self.sw_HR:str = line[1]
                    else:
                        self.sw_HR:str = line[0].replace("sw_HR=","")
                elif ( sw_ph and line[0][0:9] == "sw_phrun=" ):
                    if ( len(line) == 2 ):
                        self.sw_phrun:bool = bool(int(line[1]))
                    else:
                        self.sw_phrun:bool = bool(int(line[0].replace("sw_phrun=","")))
                elif ( sw_ph and line[0][0:10] == "sw_plt_ph=" ):
                    if ( len(line) == 2 ):
                        self.sw_plt_ph:str = line[1]
                    else:
                        self.sw_plt_ph:str = line[0].replace("sw_plt_ph=","")
                else:
                    pass
        else:
            pass

        args:str = sys.argv
        if ( (len(args) == 1 and self.mat == "") or "-h" in args ):
            helpmessage: str = """ *** USAGE: python lumin.py {mat} [-h] [-nc {nc}] [-exe {exe}] [-p] [-emod] [-ccd] [-ph]
                                -h : print helpmessage
                                -nc {nc} : number of cpu {20}
                                -exe {exe} : execute path {$HOME/qe/bin}
                                -p : sw_plt == True for all script {False}
                        [-emod] -pemod : sw_plt_emod == True {False}
                                -brav {brav} : Bravias type (cub, tet, hex, mono, tri, ...) {cub}
                                -dr {dratio} : distortion ratio {0.02}
                                -ndem {ndiv_emod} : number of devision for elastic moduli calc. {15}
                                -Bmod : calculate Bmod or not {False}
                                -egap : extract band gap energy {False}
                        [-ccd]  -pccd : sw_plt_ccd == True {False}
                                -Eabs {Eabs0} : absorption energy by first-principles calc. {2.146}
                                -Eem {Eem0} : emission energy by first-principles calc. {1.702}
                                -FC {EFCg} : Frank-Condon params. for ground state energy curve {0.214}
                                -dQ {deltaQ} : difference of the configuration coordinate b/w exited & ground {1.0}
                                -I0 {I0} : Intensity of the spectrum {1.0}
                                -gamma {gamma} : broadness params. of Lorentzian {1.0e-5}
                                -eg : generate intermediate structures b/w excited & ground or not {False}
                                -stateg {stateg} : name of the ground state {4A2g}
                                -statee {statee} : name of the excited state {4T2g}
                                -unit {sw_unit} : select unit of quantities (eV, cm^-1, nm) {eV}
                                -emin_ccd {emin_ccd} : energy min. for plotting {1.0}
                                -emax_ccd {emax_ccd} : energy max. for plotting {3.0}
                                -tempmin {tempmin} : temperature min. {1.0}
                                -tempmax {tempmax} : temperature max. {1.0e3}
                                -nmax {nmax} : cut-off number of summation, larger value is plausible {30}
                                -nde {ndiv_e} : number of devision for energy {1.0e6}
                                -ndt {ndiv_temp} : number of devision for temperature {1.0e3}
                                -ndeg {ndiv_eg} : number of devision for intermediate structures {12}
                        [-ph]   -pph : sw_plt_ph == True {False}
                                -ndim {ndim[1:3]} : number of dimension for generating supercell {1 1 1}
                                -HR {sw_HR} : switch for means to evaluate q_k (none, pos, force, both) {none}
                                -phrun : perform phonon calc. or not {False}
                                -nde_ph {ndiv_ph} : number of division for phonon spectrum {1.0e6}
                                -emin_ph {emin_ph} : energy min. for phonon spectrum {0.0}
                                -emax_ph {emax_ph} : energy max. for phonon spectrum {1.0}
                                -sigma {sigma} : standard derivation of Gaussian {6.0e-3} """
            print(helpmessage)
            sys.exit()

        if ( ( self.mat == "" ) and ( not args[1][0] == "-" ) ):
            self.mat: str = args[1]
        if ( "-emod" in args ):
            self.sw_run_emod:bool = True
        if ( "-ccd" in args ):
            self.sw_run_ccd:bool = True
        if ( "-ph" in args ):
            self.sw_run_ph:bool = True
        if ( "-nc" in args ):
            self.nc:int = int(args[args.index("-nc")+1])
        if ( "-exe" in args ):
            self.exe:str = args[args.index("-exe")+1]
        if ( "-p" in args ):
            self.sw_plt_emod:bool = True
            self.sw_plt_ccd:bool = True
            self.sw_plt_ph:bool = True
        if ( "-pemod" in args ):
            self.sw_plt_emod:bool = True
        if ( "-pccd" in args ):
            self.sw_plt_ccd:bool = True
        if ( "-pph" in args ):
            self.sw_plt_ph:bool = True
        """ emod """
        if ( "-brav" in args ):
            self.brav:str = args[args.index("-brav")+1]
        if ( "-dr" in args ):
            self.dratio:float = float(args[args.index("-dr")+1])
        if ( "-ndem" in args ):
            self.ndiv_emod:int = int(args[args.index("-ndem")+1])
        if ( "-Bmod" in args ):
            self.sw_Bmod:bool = True
        if ( "-egap" in args ):
            self.sw_egap:bool = True
        """ ccd """
        if ( "-Eabs" in args ):
            self.Eabs0:float = float(args[args.index("-Eabs")+1])
        if ( "-Eem" in args ):
            self.Eem0:float = float(args[args.index("-Eem")+1])
        if ( "-FC" in args ):
            self.EFCg:float = float(args[args.index("-FC")+1])
        if ( "-dQ" in args ):
            self.dQ:float = float(args[args.index("-dQ")+1])
        if ( "-I0" in args ):
            self.I0:float = float(args[args.index("-I0")+1])
        if ( "-gamma" in args ):
            self.gamma:float = float(args[args.index("-gamma")+1])
        if ( "-eg" in args ):
            self.sw_eg:bool = True
        if ( "-stateg" in args ):
            self.stateg:str = args[args.index("-stateg")+1]
        if ( "-statee" in args ):
            self.stateg:str = args[args.index("-statee")+1]
        if ( "-unit" in args ):
            self.sw_unit:str = args[args.index("-unit")+1]
        if ( "-emin_ccd" in args ):
            self.emin_ccd:float = float(args[args.index("-emin_ccd")+1])
        if ( "-emax_ccd" in args ):
            self.emax_ccd:float = float(args[args.index("-emax_ccd")+1])
        if ( "-tempmin" in args ):
            self.tempmin:float = float(args[args.index("-tempmin")+1])
        if ( "-tempmax" in args ):
            self.tempmax:float = float(args[args.index("-tempmax")+1])
        if ( "-nmax" in args ):
            self.nmax:int = int(args[args.index("-nmax")+1])
        if ( "-nde" in args ):
            self.ndiv_e:int = int(args[args.index("-nde")+1])
        if ( "-ndt" in args ):
            self.ndiv_temp:int = int(args[args.index("-ndt")+1])
        if ( "-ndeg" in args ):
            self.ndiv_eg:int = int(args[args.index("-ndeg")+1])
        """ ph """
        if ( "-ndim" in args ):
            self.ndim:int = [int(args[args.index("-ndim")+1+i]) for i in range(3)]
            self.ndim = np.array(self.ndim)
        if ( "-HR" in args ):
            self.sw_HR:str = args[args.index("-HR")+1]
        if ( "-phrun" in args ):
            self.sw_phrun:bool = True
        if ( "-nde_ph" in args ):
            self.ndiv_ph:int = int(args[args.index("-nde_ph")+1])
        if ( "-emin_ph" in args ):
            self.emin_ph:float = float(args[args.index("-emin_ph")+1])
        if ( "-emax_ph" in args ):
            self.emax_ph:float = float(args[args.index("-emax_ph")+1])
        if ( "-sigma" in args ):
            self.sigma:float = float(args[args.index("-sigma")+1])

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
        print("*")
        if ( self.sw_run_emod ):
            print("*** elastic moduli parameters ***")
            print("* brav: ", self.brav)
            print("* sw_plt_emod: ", self.sw_plt_emod)
            print("* sw_Bmod: ", self.sw_Bmod)
            print("* sw_egap: ", self.sw_egap)
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
            print("* sw_unit: ", self.sw_unit)
            print("* stateg: ", self.stateg)
            print("* statee: ", self.statee)
            print("* emin_ccd (eV): ", self.emin_ccd)
            print("* emax_ccd (eV): ", self.emax_ccd)
            print("* tempmin (K): ", self.tempmin)
            print("* tempmax (K): ", self.tempmax)
            print("* nmax: ", self.nmax)
            print("* ndiv_e: ", self.ndiv_e)
            print("* ndiv_temp: ", self.ndiv_temp)
            print("* ndiv_eg: ", self.ndiv_eg)
            print("*")
        if ( self.sw_run_ph ):
            print("*** phonon parameters ***")
            print("* ndim: ", self.ndim)
            print("* emin_ph: ", self.emin_ph)
            print("* emax_ph: ", self.emax_ph)
            print("* sigma: ", self.sigma)
            print("* sw_HR: ", self.sw_HR)
            print("* sw_phrun: ", self.sw_phrun)
            print("*")
        print("* --- Finish check parameters --- *")
        
    ### ----------------------------------------------------------------------------- ###
    def get_POSCAR(self, fn:str):
        """ get information from POSCAR """
        
        with open(fn, "r") as f:
            data:str = f.readlines()
        alat:float = float(data[1].split()[0]) / self.Bohr
        plat:float = [[float(d) for d in data[2].split()],
                      [float(d) for d in data[3].split()],
                      [float(d) for d in data[4].split()]]
        elements:str = [d for d in data[5].split()]
        nelems:int = [int(d) for d in data[6].split()]
        natm:int = sum(nelems)
        pos:float = []
        for dat in data[8:]:
            pos.append([float(d) for d in dat.split()])
        plat = np.array(plat)
        pos = np.array(pos)
        volume:float = (alat**3.) * abs(np.dot(plat[0],np.cross(plat[1],plat[2])))
        return (alat, plat, elements, nelems, natm, pos, volume)

    ### ----------------------------------------------------------------------------- ###
    def get_FORCE(self, fn:str):
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
    def Lorentzian(self, x:float):
        """ define Lorentzian """
        
        return (1./np.pi)*(0.5*self.gamma / (x**2.+(0.5*self.gamma)**2.))

    ### ----------------------------------------------------------------------------- ###
    def Gaussian(self, x:float):
        """ define Gaussian """
        
        return (1./(self.sigma*np.sqrt(2.*np.pi))) * np.exp(-x**2./(2.*self.sigma**2.))

    ### ----------------------------------------------------------------------------- ###
    def E2lambda(self, energy:float):
        """ transfrom energy to wavelength """
        
        return const.fceV2nm/energy  # nm

    ### ----------------------------------------------------------------------------- ###
    def lambda2E(self, lam:float):
        """ transform wavelength to energy """
        
        return const.fceV2nm/lam     # eV
