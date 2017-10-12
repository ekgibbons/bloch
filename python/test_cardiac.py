from __future__ import division
from __future__ import print_function

import numpy as np
from matplotlib import pyplot as plt

import bloch_cardiac


sliceThickness = 0.008;
nValues = int(10e3)
nExcitations = 260
T1 = 0.250
T2 = 0.020
returnType = 1

z = np.linspace(-sliceThickness*4.0,sliceThickness*4.0,
                nValues)

plt.figure(1)
plt.figure(2)
for ii in np.arange(1,260,20):
    signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
                                           ii,T1,T2,0.,returnType)

    
    plt.figure(1)
    plt.plot(z,signal.real)
    plt.title("Mx")
    plt.xlabel("slice distance [cm]")

    plt.figure(2)
    plt.plot(z,signal.imag)
    plt.title("My")
    plt.xlabel("slice distance [cm]")

    signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
                                           ii,T1,T2,0.,2)

    plt.figure(3)
    plt.plot(z,signal.real)
    plt.title("Mz")
    plt.xlabel("slice distance [cm]")
    

    
# # 90 bpm
# signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
#                                        nExcitations,T1,T2,90.,0)

# plt.figure(3)
# p1, = plt.plot(signal.real)
# p2, = plt.plot(signal.imag)
# plt.title("Signal Progression, 90bpm")
# plt.xlabel("Repetition number")
# plt.legend((p1,p2),("Real component","Imag component"),loc="best")

# 60bpm
signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
                                       nExcitations,T1,T2,0.,0)

plt.figure(4)
p1, = plt.plot(signal.real)
p2, = plt.plot(signal.imag)
plt.title("Signal Progression, 60bpm")
plt.xlabel("Repetition number")
plt.legend((p1,p2),("Real component","Imag component"),loc="best")

# # 0 bpm
# signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
#                                        nExcitations,T1,T2,0.,0)

# plt.figure(5)
# p1, = plt.plot(signal.real)
# p2, = plt.plot(signal.imag)
# plt.title("Signal Progression, 0bpm")
# plt.xlabel("Repetition number")
# plt.legend((p1,p2),("Real component","Imag component"),loc="best")


plt.show()
