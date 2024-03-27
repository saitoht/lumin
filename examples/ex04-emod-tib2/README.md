# Example04: Elastic Moduli of TiB2

## Prepare input files
- LUMIN
- POSCAR0
- tib2.scf_Bmod.in
- tib2.scf_Emod.in

## Execute the program
```shell-session
lumin -emod -nc 4 -p
```
If you can use a cluster machine, run the following command.
```shell-session
qsub job.sh
```

## Results

| Method | $ C_{11} $  | $ C_{12} $ | $ C_{44} $ |  $B$  |  $G$  |  $E$   | $ \Theta_D $ |
| ----   | ----        | ----       | ----       | ----  | ----  | ----   | ----         |
| DFT    |  |  |  |  |  |  |     |
| Expt.  |  |  |  |  |  |  |     |