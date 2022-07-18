#user/bin/python3
import numpy as np
import linecache
import statsmodels.tsa.api as smt
import matplotlib.pyplot as plt
from scipy.fftpack import fft
import math, sys

###################################################################
#
#  Usage: python3 cal_SED.py dev=0/1/2/3 num=10000/20000/30000/..
#  Note: this script is used for MS trajectory
#
#####################################################################

f_name = 'trj.txt'
dev = int(sys.argv[1].split("=")[1])
num = int(sys.argv[2].split("=")[1])
n_cell = 100        #define the number of unit cell in z axis
k_in_z = 101        #define the number of k points in z axis
pi = math.pi
k0 = pi/(12*1e-10)      #定义z轴倒格矢长度
potim = 5e-15         #取样时间间隔s
tao0 = potim*num
Fs = (1/dt)*1e-12    #抽样频率THz
natoms = int(linecache.getline(f_name, 1))
atoms = int(natoms/n_cell)

print(dev)

periodic_mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.0121831, 'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.99840316, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.921595, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584, 'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Tc': 98.9072, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.905452, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Pm': 144.9, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0245, 'Ac': 227.0277, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 239.0642, 'Am': 243.0614, 'Cm': 247.0704, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.083, 'Fm': 257.0591, 'Md': 258.0984, 'No': 259.101, 'Lr': 262.1097, 'Rf': 267.1218, 'Db': 268.1257, 'Sg': 269.1286, 'Bh': 274.1436, 'Hs': 277.1519, 'Mt': 278, 'Ds': 281, 'Rg': 282, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 288, 'Lv': 292}

def read_pos(file_name):
    pos, ele = [], []
    with open(file_name) as fi:
        for line in fi:
            ll = line.strip(' ').split()
            if len(ll) == 4:
                ele.append(ll[0])
                pos.append(list(map(float, ll[1:])))
    return np.array(pos), ele[:atoms]

def divide(array, n):
    result = []
    for i in range(n):
        x = array[i::n][:,0].tolist()
        y = array[i::n][:,1].tolist()
        z = array[i::n][:,2].tolist()
        result.append(x);result.append(y);result.append(z);
    return np.array(result)    # (natom*nframs)*time

def cal_vel(array1):
    velocity = []
    for i in range(np.shape(array1)[1]-1):
        v = (array1[:,i+1]-array1[:,i])*1e-10/dt
        velocity.append(v.tolist())
    return np.array(velocity).T

print('Reading in and dividing trajectory...')
position, element = read_pos(f_name)
pos_div = divide(position,natoms)
print('Done!')
print('Calculating velocities...')
if dev == 0:
    vel = pos_div[:,:]
elif dev == 1:
    vel = cal_vel(pos_div)
elif dev == 2:
    vel1 = cal_vel(pos_div)
    vel = cal_vel(vel1)
elif dev == 3:
    vel2 = cal_vel(pos_div)
    vel1 = cal_vel(vel2)
    vel = cal_vel(vel1)
else:
    print("You should give the variable 'dev' as 'python3 script.py dev=0,1,2,3'")
    exit()
print('Done!')
print('Calculating phonon spectral energy density (SED) function')
sed_f = []
for i in range(atoms):
    print(i)
    vx = vel[i*3::3*atoms,:]
    vy = vel[i*3+1::3*atoms,:]
    vz = vel[i*3+2::3*atoms,:]
    vx_fft1, vy_fft1, vz_fft1 = [], [], []
    vx_fftf, vy_fftf, vz_fftf = [], [], []
    for j in range(np.shape(vx)[0]):
        vx_fft = fft(vx[j,:].tolist())
        vy_fft = fft(vy[j,:].tolist())
        vz_fft = fft(vz[j,:].tolist())
        vx_fft1.append(vx_fft[:])
        vy_fft1.append(vy_fft[:])
        vz_fft1.append(vz_fft[:])
    for j in range(np.shape(vx_fft1)[1]):
        vx_fft2 = fft(np.array(vx_fft1)[:,j].tolist())
        vy_fft2 = fft(np.array(vy_fft1)[:,j].tolist())
        vz_fft2 = fft(np.array(vz_fft1)[:,j].tolist())
        vx_fftf.append(vx_fft2[:])
        vy_fftf.append(vy_fft2[:])
        vz_fftf.append(vz_fft2[:])
    vx_afft = np.array(vx_fftf).T/n_cell
    vy_afft = np.array(vy_fftf).T/n_cell
    vz_afft = np.array(vz_fftf).T/n_cell
    sed = np.sqrt(np.abs(vx_afft*(vx_afft.conjugate())+vy_afft*(vy_afft.conjugate())+vz_afft*(vz_afft.conjugate())))
    print(np.shape(sed))
    m = periodic_mass[element[i]]/(6.02e26)
    sed_f.append((np.array(sed[:])*m).tolist())
sed_final = np.sum(sed_f[:],axis=0)/2/tao0

n_data = np.shape(sed_final)[1]
f_data = np.shape(sed_final)[1]
f = [Fs*i/(n_data) for i in range(int(f_data))]
k_points = [k/n_cell for k in range(n_cell)]
print('Done!')
print('Writing in file')
with open("sed_k.out", "w") as wfe:
    wfe.write('k/fre'+' ')
    for item in f[:]:
        wfe.write(str(item)+' ')
    wfe.write('\n')
    for i in range(n_cell):
        wfe.write(str(k_points[i])+' ')
        for j in range(f_data):
            wfe.write(str(sed_final[i][j])+' ')
        wfe.write('\n')

print('Ploting figure')
a_sed_final = np.array(sed_final)
print(np.shape(a_sed_final))
print(np.shape(np.array(f)))
print(np.shape(np.array(k_points)))
levels = np.linspace(a_sed_final.min(), a_sed_final.max(), 7)
fig, ax = plt.subplots()
ax.contourf(np.array(k_points), np.array(f), a_sed_final.T, levels=levels)
plt.show()
plt.savefig('SED.png')
print('All done!!!')