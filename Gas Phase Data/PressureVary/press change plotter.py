
import numpy as np
import matplotlib.pyplot as plt

Temp,Press,H2,C2H2,C2H6,NH2,NH3,HCN,CO,H2O,H2CO,CH3OH,NO=np.genfromtxt("DiffPressure150K.csv",delimiter=",",unpack=True)


plt.plot(Press,NH2,label="NH2")
plt.xlabel("Pressure (kPa)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.savefig("NH2PressVary150K.png")
plt.show()
plt.cla()


plt.plot(Press,HCN,label="HCN")
plt.xlabel("Pressure (kPa)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.savefig("HCNPressVary150K.png")
plt.show()
plt.cla()

plt.plot(Press,NH3,label="NH3")
plt.xlabel("Pressure (kPa)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.savefig("NH3PressVary150K.png")
plt.show()
plt.cla()

plt.plot(Press,CO,label="CO")
plt.xlabel("Pressure (kPa)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.savefig("COPressVary150K.png")
plt.show()
plt.cla()
