import numpy as np
import matplotlib.pyplot as plt

import bloch_fse
import varflip
import mydisplay

mydisplay.SetPlot()

nValues = 3000
etl = 70
initialPhase = 0 # np.pi/2
returnType = 1
T1 = 1300.e3
T2 = 35.e3
esp = 4.1e3
espLRX = 2*3.6e3
nMin = 6
nLope = 13
flipFirst = 130.
flipMin = 55.
flipLope = 60.
flipLast = 100.
flipMax = 140.
dbgFlag = 0

flipsVariable = varflip.GenerateFlips(etl,esp,nMin,nLope,
                                      flipFirst,flipMin,flipLope,flipLast,flipMax,
                                      dbgFlag)


phases = np.linspace(0,np.pi/2,3)

flipsAlsop = 60.*np.ones((etl))
flipsAlsop[0] = 142.2
flipsAlsop[1] = 94.9
flipsAlsop[2] = 69.2
flipsAlsop[3] = 63.0
flipsAlsop[4] = 60.2


counter = 1
plt.figure(figsize=(20,6))
for initialPhase in phases:
    s1 = bloch_fse.BlochFSEAlsop(flipsAlsop.astype(float), nValues,
                                 etl, initialPhase, 0, T1, T2, esp)
    s2 = bloch_fse.BlochFSEGibbons(flipsVariable.astype(float), nValues,
                                   etl, initialPhase, 0, T1, T2, esp)
    # s3 = bloch_fse.BlochFSELRX(nValues, etl, initialPhase, 0, T1, T2, espLRX)
    
    plt.subplot(1,3,counter)        
    p1, = plt.plot(abs(s1))
    p2, = plt.plot(abs(s2))
    # p3, = plt.plot(abs(s3))

    if counter is 1:
        plt.ylabel("signal")
    if counter is 2:
        plt.xlabel("echo number")
    plt.title("Starting phase:   ",fontsize=24)
    plt.grid()
    plt.xlim(0, etl)
    plt.ylim(0, 0.08)
    plt.legend((p1,p2),("Alsop","Gibbons"))
    counter += 1 
plt.savefig("muscle_simulation.pdf", bbox_inches='tight')           


    
plt.show()



