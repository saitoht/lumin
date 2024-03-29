# Example03: 1D-CCD model of Ruby
A ruby, or Cr$^{3+}$-doped $\mathrm{\alpha}$-Al$_2$O$_3$,  is a well-known red gemstone.
This material is used for the first solid-state laser by Maiman in 1960.

We simulate the optical spectra of this system using first-principles calculations. 

## Prepare input files
- LUMIN
- POSCAR_4A2g
- POSCAR_4T2g

## Execute the program
For optical spectra in 1D-CCD model.
```shell-session
lumin -ccd -p
```

## Results
Output text about spectra
- lumin.results_spec

Figures
- Spectrum.pdf
- Epeak_Temp.pdf
- Stokes_Temp.pdf
- FWHM_Temp.pdf