# ChipSum.MD
ChipSum's toolkits serve computational chemistry applications in CDCS, including molecular dynamics simulation and the calculation of the physical and chemical properities.

   - **Preproc:** Some scripts to deal with the problems occur before calculation or simulation, such as the input parameters suggestion, format transformation and so on.
   - **Postproc:** Some scripts to deal with the problems occur after calculation or simulation, such as the data analyse, plotting graph and so on.

## Functions

### Mean Squared Displacement (MSD)
`MSDVASP.py` is used to calculate MSD from VASP molecular dynamic simulation output file `XDATCAR`.
`MSDQE.py` is used to calculate MSD from QE molecular dynamic simulation output file `.pos`.<br>

   - **Usage:** `python3 MSDVASP.py/MSDQE.py`
   - **Note:** `XDATCAR` file shoule be in the directory. The detial notes are in the scripts. 

### Bond length
`bondVASP.py` is used to calculate and analyze the bond length from VASP molecular dynamic simulation output file `XDATCAR`.<br>

   - **Usage:** `python3 bondVASP.py ele1 ele2`
   - **Note:** `ele1` and `ele2` are element symblo case-insensitive. The detial notes are in the scripts.

### Fast Fourier Transform (FFT)
`FFT.py` is used to make FFT of atom position from `MSDVASP.py` and `MSDQE.py` output file `.pos`.<br>

### Radial Distribution Function (RDF)
`RDFVASP.py` is used to calculate RDF of pair atoms from `MSDVASP.py` and `MSDQE.py` output file `.pos`.<br>

   - **Usage:** `python3 RDFVASP.py`

### Plot motion trajectory graphs
`plot_trajectory.py` is used to plot motion trajectory graphs of atoms from MSDVASP.py and MSDQE.py output file .pos.<br>

### Spectral Energy Density (SED)
`cal_SED.py` is used to calculate the SED from MD trajectory and then Phonon dispersion and Phonon density of states (DOS).<br>

   - **Note:** Details notes is in the directory of `Postproc/SED/`
   - **Example:** SWNT

![SWNT](https://github.com/EltonYH/ChipSum.MD/blob/main/Postproc/img/swnt_small.png)<br>

### Velocity Autocorrelation Function (VACF)
The VACF is an important tool to analyze the atom vibrations in the MD simulations. The script `cal_VACF.py` will output the VACF of different atoms seperated by atom ID and three axis. Furthermore, it will also calculated the FFT of VACF and give the Lorentz distribution fitting of it. Besides, the force constants of all the atoms in your system will be also calculated.

   - **Note:** Details notes is in the directory of `Postproc/VACF/`
