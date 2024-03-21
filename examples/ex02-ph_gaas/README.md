# Example02: Phonon Spectrum of GaAs
We calculate the phonon spectrum of GaAs by using Quantum ESPRESSO.

## Prepare input files
- LUMIN
- gaas.scf.in
- header.in
- band.conf

## Execute the program
This is a bit heavy calculation, so I recommend you to use cluster machine to calculate, or reduce the computation costs.

To submit the job, you will execute the following command. 
```shell-session
qsub job.sh
```

## Results
Output text
- lumin.out

After DFT calculations, you can see the results of phonon spectrum by the following command.
```shell-session
phonopy -p band.conf
```