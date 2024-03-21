# Example01: Elastic Moduli of Diamond

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

## Results

| Method | $ C_{11} $  | $ C_{12} $ | $ C_{44} $ |  $B$  |  $G$  |  $E$   | $ \Theta_D $ |
| ----   | ----        | ----       | ----       | ----  | ----  | ----   | ----         |
| DFT    | 1057 | 122 | 576 | 434 | 530 | 1134 | 2232    |
| Expt.  | 1077 | 125 | 577 | 442 | 534 | 1146 | 1883    |