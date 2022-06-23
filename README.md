# ChipSum.MD
ChipSum's toolkits serve computational chemistry applications in CDCS, including molecular dynamics simulation and the calculation of the physical and chemical properities.

MSDVASP.py: Calculate Mean Squared Displacement from VASP molecular dynamic simulation output file XDATCAR.

MSDQE.py: Calculate Mean Squared Displacement from QE molecular dynamic simulation output file .pos.

bondVASP.py: Calculate and analyze bond length from VASP molecular dynamic simulation output file XDATCAR.

FFT.py: Make FFT(Fast Fourier Transform) of atom position from MSDVASP.py and MSDQE.py output file .pos.

RDFVASP.py: Calculate RDF(Radial distribution function) of pair atoms from MSDVASP.py and MSDQE.py output file .pos.

plot_trajectory.py: Plot motion trajectory graphs of atoms from MSDVASP.py and MSDQE.py output file .pos.

cal_SED.py: Calculate the spectral energy density (SED) from MD trajectory and then Phonon dispersion and Phonon density of states (DOS).
