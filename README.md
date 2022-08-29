# ChipSum.MD
ChipSum's toolkits serve computational chemistry applications in CDCS, including molecular dynamics simulation and the calculation of the physical and chemical properities.
        

## Function

### Mean Squared Displacement (MSD)
`MSDVASP.py` is used to calculate MSD from VASP molecular dynamic simulation output file `XDATCAR`.
`MSDQE.py` is used to calculate MSD from QE molecular dynamic simulation output file `.pos`.<br>

   - ** Usage: ** `python3 MSDVASP.py/MSDQE.py`

### Bond length
`bondVASP.py`: Calculate and analyze bond length from VASP molecular dynamic simulation output file `XDATCAR`.<br>

### Fast Fourier Transform (FFT)
`FFT.py`: Make FFT of atom position from `MSDVASP.py` and `MSDQE.py` output file `.pos`.<br>

### Radial Distribution Function (RDF)
`RDFVASP.py`: Calculate RDF of pair atoms from `MSDVASP.py` and `MSDQE.py` output file `.pos`.<br>

### Plot motion trajectory graphs
`plot_trajectory.py`: Plot motion trajectory graphs of atoms from MSDVASP.py and MSDQE.py output file .pos.<br>

### Spectral Energy Density (SED)
`cal_SED.py`: Calculate the SED from MD trajectory and then Phonon dispersion and Phonon density of states (DOS).<br>
![SWNT](https://github.com/EltonYH/ChipSum.MD/blob/main/Postproc/img/swnt_small.png)
