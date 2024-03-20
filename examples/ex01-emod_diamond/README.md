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