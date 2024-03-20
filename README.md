# lumin
Simulate luminescent properties of materials using the information from first-principles calculations

## Requirement
For first-principles calculations
- Quantum ESPRESSO

For Phonon calculations
- Phonopy

For Python modules
- numpy
- matplotlib
- pathlib
- scipy
- subprocess
- yaml
- sklearn

## Usage
Before you use the program "lumin", you should activate the path by writing
```shell-session
export PATH=$PATH:$HOME/lumin
```
in ~/.bashrc file. 
After that, you should execute
```shell-session
source ~/.bashrc
```
to reflect the settings or exit the terminal once.

In the case of ZSH shell, you will see ~/.zshrc instead of ~/.bashrc.
You can check which SHELL you are using by the following command.
```shell-session
echo $SHELL
```

For elastic moduli calculations
```shell-session
lumin -emod
```

I checked the results for cubic system.
You should be careful if you calculate non-cubic system.

For phonon calculations
```shell-session
lumin -ph
```

For configuration coordinate diagram calculations
```shell-session
lumin -ccd
```

## Tutorials
See examples directory
- Ex01: Elastic Moduli of Diamond
- Ex02: Phonon Calculation of GaAs
- Ex03: Configuration Coordinate Model of Ruby

## References
For configuration coordinate diagram
- M. A. Reshchikov & H. Morkoc, Luminescence properties of defects in GaN, J. Appl. Phys. 97 (2005) 061301.
- Y. Jia, A. Miglio, S. Ponce, M. Mikami, & X. Gonze, First-principles study of the luminescence of Eu2+-doped phosphors, PRB 96 (2017) 125132.