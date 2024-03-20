""" Lumin.py: Estimation of luminescent properties of materials
           Consolidate cls_emod.py, cls_ccd.py, and cls_ph.py """
import sys, time
import cls_subs as subs
import cls_emod as emod
import cls_ph as ph
import cls_ccd as ccd

if __name__=="__main__":
    prms = subs.subs()
    now = time.strftime("%Y/%m/%d %H:%M:%S")
    tmbf = time.time()
    print("* --- lumin.py --- *")
    print("*** JOB Start at {now}.".format(now=now))
    print("*")
    if ( prms.sw_run_emod == True ):
        """ Elastic moduli """
        emod.emod()
    if ( prms.sw_run_ph == True ):
        """ Phonon """
        ph.ph()
    if ( prms.sw_run_ccd == True ):
        """ configuration coordinate model """
        ccd.ccd()
    tmaf = time.time()
    print("* --- All calculations are finished! --- *")
    print("* Computation time: {dt} (sec)".format(dt=tmaf-tmbf))
    print("*")
    sys.exit()
