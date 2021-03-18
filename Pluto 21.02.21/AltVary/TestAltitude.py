import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline

Altitudea,SSDa,Tempa,Pressa,Rada,Metha,CO2a,H2a,CH4a,C2H2a,C2H6a,N2a,NH3a,CNa,HCNa,COa,H2Oa,H2COa,CH3OHa,NOa=np.genfromtxt("Triala.csv",delimiter=",",unpack=True)
Altitudeb,SSDb,Tempb,Pressb,Radb,Methb,CO2b,H2b,CH4b,C2H2b,C2H6b,N2b,NH3b,CNb,HCNb,COb,H2Ob,H2COb,CH3OHb,NOb=np.genfromtxt("Trialb.csv",delimiter=",",unpack=True)
Altitudec,SSDc,Tempc,Pressc,Radc,Methc,CO2c,H2c,CH4c,C2H2c,C2H6c,N2c,NH3c,CNc,HCNc,COc,H2Oc,H2COc,CH3OHc,NOc=np.genfromtxt("Trialc.csv",delimiter=",",unpack=True)
Altituded,SSDd,Tempd,Pressd,Radd,Methd,CO2d,H2d,CH4d,C2H2d,C2H6d,N2d,NH3d,CNd,HCNd,COd,H2Od,H2COd,CH3OHd,NOd=np.genfromtxt("Triald.csv",delimiter=",",unpack=True)
sort=Altitudeb.argsort()
Tempb=Tempb[sort]
Altitudeb=Altitudeb[sort]
Pressb=Pressb[sort]

xnew=np.linspace(Altitudeb.min(),Altitudeb.max(),100)
spl=make_interp_spline(Altitudeb,Tempb,k=3)
smooth=spl(xnew)



fig,ax1=plt.subplots()
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Altitude (km)")
ax1.plot(Tempb, Altitudeb)
ax2=ax1.twinx()
ax2.plot(Tempb,Pressb,alpha=0)
ax2.set_ylabel("Pressure (kPa)")
ax2.set_yscale("log")
ax2.invert_yaxis()
plt.savefig("AltitudeTempPressValuesPlu.png",bbox_inches="tight")
plt.show()
plt.cla()

fig,ax1=plt.subplots()
ax1.set_ylabel("Temperature (K)")
ax1.set_xlabel("Altitude (km)")
ax1.plot(xnew,smooth)
ax2=ax1.twinx()
ax2.plot(Altitudeb,Pressb,alpha=0)
ax2.set_ylabel("Pressure (kPa)")
ax2.set_yscale("log")
ax2.invert_yaxis()
plt.show()
plt.cla()