### Lumin.py: Estimation of luminescent properties of materials
###           Consolidate cls_emod.py, cls_ccd.py, and cls_ph.py
import sys
import cls_subs as subs
import cls_emod as emod
import cls_ph as ph
import cls_ccd as ccd

if __name__=="__main__":
    const = subs.subs()
    prms = subs.get_prms()
    print("* --- lumin.py --- *")
    print("*")
    if ( prms.sw_run_emod == True ):
        # Elastic moduli
        emod.run_emod()
    if ( prms.sw_run_ph == True ):
        # Phonon
        ph.run_phonon()
    if ( prms.sw_run_ccd == True ):
        # configuration coordinate model
        ccd.run_ccd()
    print("* --- All calculations are finished! --- *")
    print("*")
    sys.exit()
