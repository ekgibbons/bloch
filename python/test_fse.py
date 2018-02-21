import numpy as np
import matplotlib.pyplot as plt

from mri import blochsequence
# import varflip
# import mydisplay

# mydisplay.SetPlot()

number = 3

nValues = 5000
etl = 76
initialPhase = 0 # np.pi/2
returnType = 1
T1 = 1300.e3
T2 = 35.e3
esp = 4.1e3
espLRX = 2*3.6e3
nMin = 6
nLope = 11
flipFirst = 145.
flipMin = 55.
flipLope = 60.
flipLast = 100.
flipMax = 140.
dbgFlag = 0

flipsAlsop = 60.*np.ones((etl))
flipsAlsop[0] = 142.2
flipsAlsop[1] = 94.9
flipsAlsop[2] = 69.2
flipsAlsop[3] = 63.0
flipsAlsop[4] = 60.2

flipsVariable = flipsAlsop
# flipsVariable = varflip.GenerateFlips(etl,esp,nMin,nLope,
#                                       flipFirst,flipMin,flipLope,flipLast,flipMax,
#                                       dbgFlag)



etl = 0
initialPhase = np.pi/4

s1 = blochsequence.FSEAlsop(flipsAlsop.astype(float), nValues,
                            etl, initialPhase, T1, T2, esp, 1)
s1z = blochsequence.FSEAlsop(flipsAlsop.astype(float), nValues,
                             etl, initialPhase, T1, T2, esp, 2)


s2 = blochsequence.FSEGibbons(flipsVariable.astype(float), nValues,
                              etl, initialPhase, T1, T2, esp, 1)
s2z = blochsequence.FSEGibbons(flipsVariable.astype(float), nValues,
                               etl, initialPhase, T1, T2, esp, 2)


z = np.linspace(-3*5,3*5,nValues)

plt.figure(figsize=(22,6))
plt.subplot(1,3,1)
plt.plot(z,s1.real)
# plt.plot(z[2000:3000],s1[2000:3000].real)
plt.xlabel("z-location [mm]")
#plt.title("Mx")
plt.grid()
plt.ylim(-1.1,1.1)

plt.subplot(1,3,2)
plt.plot(z,s1.imag)
# plt.plot(z[2000:3000],s1[2000:3000].imag)
plt.xlabel("z-location [mm]")
#plt.title("My")
plt.grid()
plt.ylim(-1.1,1.1)

plt.subplot(1,3,3)
plt.plot(z,s1z.real)
# plt.plot(z[2000:3000],s1z[2000:3000].real)
plt.xlabel("z-location [mm]")
#plt.title("My")
plt.grid()
plt.ylim(-1.1,1.1)
plt.savefig("alsop_%i.pdf" % number,bbox_inches="tight")

print("alsop net magnetization: %f + j%f" % (np.mean(s1).real,np.mean(s1).imag))

plt.figure(figsize=(22,6))
plt.subplot(1,3,1)
plt.plot(z,s2.real)
# plt.plot(z[2000:3000],s2[2000:3000].real)
plt.xlabel("z-location [mm]")
#plt.title("Mx")
plt.grid()
plt.ylim(-1.1,1.1)

plt.subplot(1,3,2)
plt.plot(z,s2.imag)
# plt.plot(z[2000:3000],s2[2000:3000].imag)
plt.xlabel("z-location [mm]")
#plt.title("My")
plt.grid()
plt.ylim(-1.1,1.1)

plt.subplot(1,3,3)
plt.plot(z,s2z.real)
# plt.plot(z[2000:3000],s2z[2000:3000].real)
plt.xlabel("z-location [mm]")
#plt.title("Mz")
plt.grid()
plt.ylim(-1.1,1.1)
plt.savefig("gibbons_%i.pdf" % number,bbox_inches="tight")

print("gibbons net magnetization: %f + j%f" % (np.mean(s2).real,np.mean(s2).imag))


# counter = 1
# plt.figure(figsize=(16,12))
# for initialPhase in phases:
#     s1 = blochsequence.BlochFSEAlsop(flipsAlsop.astype(float), nValues,
#                                  etl, initialPhase, 0, T1, T2, esp)
#     s2 = blochsequence.BlochFSEGibbons(flipsVariable.astype(float), nValues,
#                                    etl, initialPhase, 0, T1, T2, esp)
#     s3 = blochsequence.BlochFSELRX(nValues, etl, initialPhase, 0, T1, T2, espLRX)
    
#     plt.subplot(2,2,counter)        
#     p1, = plt.plot(abs(s1))
#     p2, = plt.plot(abs(s2))
#     p3, = plt.plot(abs(s3))
#     if counter is 3 or counter is 4:
#         plt.xlabel("echo number")
#     plt.ylabel("signal")
#     plt.title("Starting phase: %.003f" % initialPhase)
#     plt.grid()
#     plt.xlim(0, etl)
#     plt.ylim(0, 0.08)
#     plt.legend((p1,p2,p3),("Alsop","Gibbons","LRX"))
#     counter += 1 
# plt.savefig("muslce_simulation.pdf", bbox_inches='tight')           

# sGibbons = blochsequence.BlochFSEGibbons(flipsVariable.astype(float), nValues,
#                                      etl, np.pi/4, 2, T1, T2, esp)

# z = np.linspace(-3*5,3*5,nValues)
# plt.figure()
# plt.plot(z,sGibbons.real)
# plt.xlabel("slice position [mm]")
# plt.ylabel("Mz")
# plt.ylim(-0.5, 1.1)
# plt.title("Gibbons")
    
# sCPMG = blochsequence.BlochFSECPMG(flipsVariable.astype(float), nValues, etl, 0, T1, T2, esp)

# plt.figure()
# plt.plot(abs(sCPMG))
# plt.title("CPMG")

    
plt.show()







# plt.plot(s1.real)
# plt.plot(s1.imag)
# plt.plot(s2.real)

# sCPMG = blochsequence.BlochFSECPMG(flipsVariable.astype(float), nValues, etl, 0, T1, T2, esp)

# plt.plot(abs(sCPMG))
    
# s2 = blochsequence.BlochFSEAlsop(flipsVariable.astype(float), nValues, 1, initialPhase, 1, T1, T2, esp)
# s3 = blochsequence.BlochFSEAlsop(flipsVariable.astype(float), nValues, 1, initialPhase, 2, T1, T2, esp)


# plt.figure()
# plt.plot(s2.real)
# plt.plot(s2.imag)
# plt.plot(s3.real)
plt.show()
