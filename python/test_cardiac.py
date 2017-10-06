from __future__ import division
from __future__ import print_function

import numpy as np
from matplotlib import pyplot as plt

import bloch_cardiac

sliceThickness = 8./1000;
nValues = int(10e3)
nExcitations = 260
T1 = 1.0
T2 = 0.035
returnType = 0

signal = bloch_cardiac.BlochCardiacSat(sliceThickness,nValues,
                                       nExcitations,T1,T2,returnType)

plt.figure()
plt.plot(signal.real)
plt.plot(signal.imag)
plt.show()
