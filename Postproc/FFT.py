#/usr/bin/python3
import numpy
import math,cmath
import matplotlib.pyplot as plt
###################################################################################
#              Fast Fourier Transformation                                        #
#              N                                                                  #
#      X(k) = sum x(n)*exp(-i*2pi*(k-1)*(n-1)/N)), 1<=k<=N                        #
#             n=1                                                                 #
#                                                                                 #
#  1. 'N' is the data number you will transform.                                  #
#  2. 'dt' is the time interval between each position, unit s.                    #
#  3. You must modify the input .pos file in line 27.                             #
#Note:                                                                            #
#  1. This code is now only approprate for atomic type whose atomic number is 1.  #
#  2. This code is suggested to use after MSDVASP.py/MSDQE.py.                    #
###################################################################################

N = 5000                               #样品个数
dt = 1.0E-15*1.0                     #取样时间间隔(s)

pi = 3.1415926536
Fs = 1/dt*1E-12                        #取样频率(THz)
DFT_X, DFT_Y, DFT_Z = [], [], []              #存储傅里叶级数
pos = []                                #坐标信息(A)
f = []
with open('Fe.pos') as file:
  for line in file:
    el = list(line.strip().split(' '))
    po =[]
    for element in el[1:]:
      po.append(float(element))
    pos.append(po)
cpos = numpy.matrix(pos);
if N <= len(pos):
  mpos = cpos[-N:,:]
elif N > len(pos):
  zeros = numpy.zeros([N-len(pos),3],dtype=float)
  mpos = numpy.vstack((cpos,zeros))

ed = numpy.zeros([N,N],dtype=complex)
print(numpy.shape(mpos))
for k in range(0,N):
  if k%100 == 0: print(k)
  for l in range(k,N):
    t = numpy.exp(-1j*2*pi*k*l/N);
    ed[l,k] = t; ed[k,l] = t;
  f.append(Fs*k/N)

med = numpy.matrix(ed)
DFT = med*mpos/N*2
aDFT = numpy.array(abs(DFT))
print(numpy.shape(aDFT))

file_name = 'DFT.out'
with open(file_name,"w") as file_object:
  file_object.write("f/THz   pDFT_X   pDFT_Y   pDFT_Z\n")
  for i in range(0,int(numpy.shape(aDFT)[0]/2)+1):
    file_object.write(str(f[i])+" ")
    for j in range(0,numpy.shape(aDFT)[1]):
      file_object.write(str(aDFT[i,j])+" ")
    file_object.write("\n")
