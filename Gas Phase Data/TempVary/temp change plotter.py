# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

Temp,Press,H2,C2H2,C2H6,NH2,NH3,HCN,CO,H2O,H2CO,CH3OH,NO=np.genfromtxt("DiffTemp300kPa.csv",delimiter=",",unpack=True)



plt.plot(Temp,NH2,label="NH2")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.savefig("NH2TempVary300kPA.png")
plt.show()
plt.cla()

plt.plot(Temp,HCN,label="HCN")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.savefig("HCNTempVary300kPa.png")
plt.show()
plt.cla()

plt.plot(Temp,NH3,label="NH3")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.savefig("NH3TempVary300kPa.png")
plt.show()
plt.cla()

plt.plot(Temp,CO,label="CO")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Concentration, n (cm$^{-3}$)")
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.savefig("COTempVary300kPa.png")
plt.show()
plt.cla()