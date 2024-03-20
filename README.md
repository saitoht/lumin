# lumin
Simulate luminescent properties of materials using the information from first-principles calculations

## Requirement
For first-principles calculations
- Quantum ESPRESSO

For Python modules
- numpy
- matplotlib
- pathlib
- scipy
- subprocess
- yaml
- sklearn

## Usage
For elastic moduli calculations
```shell-session
lumin -emod
```

For phonon calculations
```shell-session
lumin -ph
```

for configuration coordinate diagram calculations
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