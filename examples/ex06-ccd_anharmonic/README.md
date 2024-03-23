# Example06: Anharmonic effects in 1D-CCD model of Ruby
We simulate the optical spectra of this system using ${\Delta}$SCF method.

## Prepare input files
- LUMIN
- POSCAR_4A2g
- POSCAR_4T2g
- header_4A2g.in
- header_4T2g.in

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
- lumin.results_spec

Figures
