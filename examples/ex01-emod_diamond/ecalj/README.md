# Example01: Elastic Moduli of Diamond (Primitive Cell)

## Prepare input files
- LUMIN
- POSCAR0
- c_dia.scf_Bmod.in
- c_dia.scf_Emod.in
- C.pbe-n-kjpaw_psl.1.0.0.UPF

## Execute the program
```shell-session
lumin -emod -nc 4 -p
```
If you can use a cluster machine, run the following command.
```shell-session
qsub job.sh
```

In order to see the results of fitting
```shell-session
lumin -emod -p
```

## Results
Output file
- lumin.out

| Method | ${ C_{11} }$ (GPa)  | ${ C_{12} }$ (GPa) | ${ C_{44} }$ (GPa) |  ${B}$ (GPa)  |  ${G}$ (GPa)  |  ${E}$ (GPa)  | ${ \Theta_D }$ (K) |
| ----   | ----        | ----       | ----       | ----  | ----  | ----   | ----         |
| DFT    | 1072 | 136 | 576 | 448 | 530 | 1145 | 2234    |
| Expt.  | 1077 | 125 | 577 | 442 | 534 | 1146 | 1883    |

Figures
- Bmod_fit.pdf
- cub-uni_fit.pdf
- cub-xy_fit.pdf
- cub-vol_fit.pdf