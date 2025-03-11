# Example06: Anharmonic effects in 1D-CCD model of Ruby
We simulate the optical spectra of this system using ${\Delta}$SCF method.

## Prepare input files
- LUMIN
- POSCAR_4A2g
- POSCAR_4T2g
- header_4A2g.in
- header_4T2g.in
- pot/Al.pbesol-n-kjpaw_psl.1.0.0.UPF
- pot/Cr.pbesol-spn-kjpaw_psl.1.0.0.UPF
- pot/O.pbesol-n-kjpaw_psl.0.1.UPF

## Execute the program
Here, in order to calculate the intermediate states automatically, you will activate "sw_eg" in LUMIN as
```shell-session
 sw_eg=1
```
, and run the following command.
```shell-session
lumin -ccd
```

I strongly recommend to use cluster machine for this example
```shell-session
qsub job.sh
```

For optical spectra in 1D-CCD model.
```shell-session
lumin -ccd -p
```

## Results
Output text about spectra
- lumin.out
- DATA_4A2g
- DATA_4T2g

Figures
- figs/1DCCD_4T2g.pdf
- figs/Ecurve_fit_DATA_4T2g.pdf
- figs/Epeak_Temp.pdf
- figs/FWHM_Temp.pdf
- figs/Spectrum.pdf
- figs/Stokes_Temp.pdf
