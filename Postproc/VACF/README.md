# Velocity Autocorrelation Function (VACF)
The VACF is an important tool to analyze the atom vibrations in the MD simulations. The script `cal_VACF.py` will output the VACF of different atoms seperated by atom ID and three axis. Furthermore, it will also calculated the FFT of VACF and give the Lorentz distribution fitting of it. Besides, the force constants of all the atoms in your system will be also calculated.

## Usage 
```shell
python3 cal_VACF.py dev=0/1/2/3 num=10000/20000/30000/...
```
`dev`: the derivation order of trajectory position.<br>
`num`: the number of MD last steps used to calculate VACF.<br>

## Note 
1. It is suggested to use this script after `output_pos.py` and `comb.sh`
2. The script `output_pos.py` is to seperate the XDATCAR by different atoms and the positions are converted to absolute coordinates.
3. The script `comb.sh` is to combine all the seperated `.pos` files.

