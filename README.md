# lumin
Simulate luminescent properties of materials using the information from first-principles calculations.

## Requirement
For first-principles calculations
- Quantum ESPRESSO
- (ecalj)

For Phonon calculations
- Phonopy
- (Quantum ESPRESSO)
- (alamode)

For Python modules
- numpy
- matplotlib
- pathlib
- scipy
- subprocess
- yaml
- sklearn

## Installation
Type the following command to install.
```shell-session
git clone https://github.com/saitoht/lumin.git
```

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
Then, "lumin" command will be activated.

In the case of ZSH, you will see ~/.zshrc instead of ~/.bashrc.
You can check which SHELL you are using by the following command.
```shell-session
echo $SHELL
```

## Usage
For elastic moduli calculations
```shell-session
lumin -emod
```

I checked the results for cubic & hexagonal system.
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
- Ex04: Elastic Moduli of TiB2
- (Ex05: Multi-D CCD model from Phonon spectrum)
- Ex06: Anharmonic effects in 1D-CCD model
- (Ex07: Automatic evaluation of Elastic Moduli & Energy Band Gap)
- (Ex08: Raman process?)
- (Ex09: Detailed Analysis on Adiabatic Energy Curve)
- (Ex10: Luminescence from Charged Defects)
- (Ex11: Nephelauxetic effects)

## References
For elastic moduli calculations
- M. Jamal, S. J. Asadabadi, I. Ahmad, H. A. R. Aliabad, Elastic constants of cubic crystals, Computational Materials Science 95 (2014) 592-599.
- A. H. Reshak, M. Jamal, DFT Calculation for Elastic Constants of Tetragonal Strucrure of Crystalline Solids with WIEN2k Code: A New Package (Tetra-elastic), Int. J. Electrochem. Sci., 8 (2013) 12252.
- Z. Zhang, Z. H. Fu, R. F. Zhang, D. Legut, and H. B. Guo, Anomalous mechanical strengths and shear deformation paths of Al2O3 polymorphs with high ionicity, RCS Advances.
- K. B. Panda, K. S. Ravi Chamdram, Determination of elastic constants of titanium diboride from first principles using FLAPW implementation of the density functional theory, Comp. Mater. Sci. (2006) 134.
- Y. Zhou, A. M. Tehrani, A. O. Oliynyk, A. C. Duke, & J. Brogch, Identifying an efficient, thermally robust inorganic phosphor host via machine learning, Nature Commn. 9 (2018) 4377.

For phonon calculations
- A. Alkauskas, B. B. Buckley, D. D. Awschalom, C. G. Van de Walle, First-principles theory of the luminescence for the triplet transition in diamond NV centres, New J. Phys. 16 (2014) 073026.

For configuration coordinate diagram
- M. A. Reshchikov & H. Morkoc, Luminescence properties of defects in GaN, J. Appl. Phys. 97 (2005) 061301.
- Y. Jia, A. Miglio, S. Ponce, M. Mikami, & X. Gonze, First-principles study of the luminescence of Eu2+-doped phosphors, PRB 96 (2017) 125132.

## Things to do
- ecalj version for emod
- automatic calculation of elastic modulus & band gap energy
- optical spectra at specific temperature, plot figures (correspondence in color & wavelength)
- Excitons, Polarons & Polaritons???
- Jahn-Teller effects?
- connection to eigloc code
