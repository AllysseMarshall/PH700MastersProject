
import numpy as np
import matplotlib.pyplot as plt

CH4,H2,C2H2,C2H6,NH2,NH3,CN,HCN,CH3CN,C6H6,H2CO,CH3OH,NO=np.genfromtxt("CH4SurfFreezeCH3CNNormalised.csv",delimiter=",",unpack=True)

 
plt.plot(CH4,(HCN/np.median(CH3CN))*9000,label="HCN")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeHCNnormalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(NH2/np.median(CH3CN))*9000,label="NH2")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeNH2normalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(NH3/np.median(CH3CN))*9000,label="NH3")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeNH3normalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(CN/np.median(CH3CN))*9000,label="CN")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeCNnormalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(CH3CN/np.median(CH3CN))*9000,label="CH3CN")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeCH3CNnormalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(C6H6/np.median(CH3CN))*9000,label="C6H6")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u))")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeC6H6normalmed.png")
plt.show()
plt.cla()

plt.plot(CH4,(C2H2/np.median(CH3CN))*9000,label="C2H2")
plt.xlabel("CH4 Composition (%)")
plt.ylabel(r"Intensity, (a.u)")
plt.legend()
#plt.yscale("log")
plt.savefig("CH4ChangeC2H2normalmed.png")
plt.show()
plt.cla()