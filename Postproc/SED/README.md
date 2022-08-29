# Spectral Energy Density (SED) 
SED is the phonon energy spectral density, which is used to analyze the phonon dispersion of materials after MD simulations.<br>

## Usage
`python3 cal_SED.py dev=0/1/2/3 num=10000/20000/30000/..`<br>

`dev`: the order of derivation of atom position, 1st order means velocity, 2nd order means acceleration....<br>
`num`: the number of steps you will used to calculate the SED, which means the last num MD steps will be adopted. <br>

Note: this script is used for MS trajectory, thus, you should convert your MS trajectory file (xtd) to reading file (xyz in txt).<br>
You can use the script of "xtd2xyz.pl" which will be found in the internet, and the author is Andrea Minoia.