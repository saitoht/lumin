# Example01: 1D-CCD model of Ruby

## Prepare input files
- LUMIN
- POSCAR_4A2g
- POSCAR_4T2g

## Execute the program
For optical spectra in 1D-CCD model.
```shell-session
lumin -ccd -p
```

If you want to calculate the intermediate states automatically, you will activate "sw_eg" in LUMIN as
```shell-session
 sw_eg=1
```
, and run the following command.
```shell-session
lumin -ccd
```

## Results
