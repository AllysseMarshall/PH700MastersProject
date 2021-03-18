import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

def PositionFunc(species, mass):
    print(np.where(species[2] == mass))



species=np.genfromtxt("species list.csv",delimiter=",",unpack=True)

#PositionFunc(species, 52.07416)
konc=np.genfromtxt("koncTit.csv",delimiter=",",unpack=True, usecols=(tuple([int(i) for i in species[0]])))
konc2=np.genfromtxt("konct2.csv",delimiter=",",unpack=True, usecols=(tuple([int(i) for i in species[0]])))
#konc3=np.genfromtxt("konc2psurf.csv",delimiter=",",unpack=True, usecols=(tuple([int(i) for i in species[0]])))
konc4=np.genfromtxt("konct2aer.csv",delimiter=",",unpack=True, usecols=(tuple([int(i) for i in species[0]])))

#nit=konc[51]
#print(nit)
#print(konc)
#print(species[0])


sort=np.argsort(species[2])
molecularmass=species[2][sort]
koncsorted=konc[sort]
#konc2sorted=konc2[sort]
#konc3sorted=konc3[sort]
#konc4sorted=konc4[sort]
molecularmassTit=molecularmass[np.where(koncsorted>1e-6)]
#molecularmass2=molecularmass[np.where(konc2sorted>1e-6)]
#molecularmass3=molecularmass[np.where(konc3sorted>1e-6)]
#molecularmass4=molecularmass[np.where(konc4sorted>1e-6)]
koncsorted=koncsorted[np.where(koncsorted>1e-6)]
#konc2sorted=konc2sorted[np.where(konc2sorted>1e-6)]
#konc3sorted=konc3sorted[np.where(konc3sorted>1e-6)]
#konc4sorted=konc4sorted[np.where(konc4sorted>1e-6)]
#print(molecularmass)


#plt.stem(molecularmass,konc\nit,label="1.5% CH4")
#plt.xlabel("Molecular Mass (amu)")
#plt.ylabel(r"Concentration, (au)")
#plt.legend()
#plt.yscale("log")
#plt.tight_layout()
#plt.savefig("TitCompAerFreNormalised.png")
#plt.show()
#plt.cla()

#plt.hist(molecularmass2,konc2,bins=len(molecularmass2),label="2% CH4")
#plt.hist2d(molecularmass3,konc3,label="3% CH4")
#plt.hist2d(molecularmass4,konc4,label= "4% CH4")
#plt.xlabel("Molecular Mass (amu)")
#plt.ylabel(r"Concentration, n (cm$^{-3}$)")
#plt.legend()
#plt.yscale("log")
#plt.tight_layout()
#plt.savefig("TitMazExpTrial.png")
#plt.show()
lims = [(0,40),(40,80),(80,105)]
norm=True
normVal = 73
labels = [[0,'H'], [1,'H2'], [3,'C'], [2,'CH'], [50,'C6H5'],[121,'HNCCHCN'], [4,'C2'],[24,'C4H'],[9,'CH4'],[80, 'N2H2'],[32,'C4H9'],[95,'CH3NH2'],[77,'NH2'],[33,'C4H10'], [46,'C6H'], [176,'CH3OOH'], [11,'C2H2'],[45,'C5H12'], [15,'C2H6'], [101,'C2N2'],[103,'C6N2'], [139,'C3H7CN'], [113,'C2H5CN'], [51,'C6H6'], [72,'N'], [73,'N2'], [78,'NH3'], [83,'CN'],[85,'HCN'],[154,'CO'],[155,'CO2'],[145,'C6H5CN'],[152,'O3']]
#[151,'O2']
for i,val in enumerate(lims):
    #if norm:
     #   markerline, stemlines, baseline = plt.stem(molecularmass3,konc3sorted/(konc3[normVal]),markerfmt='o', label= "Experimental Surfaces")
    #else:
     #   markerline, stemlines, baseline = plt.stem(molecularmass3,konc3sorted,markerfmt='o', label= "Experimental Surfaces")
    #plt.setp(stemlines, 'color', plt.getp(markerline,'color'))
    #plt.setp(stemlines, 'alpha', 0.5)
    #plt.setp(markerline, color='blue')
    #plt.setp(stemlines, color='blue')
    #plt.setp(markerline, 'alpha', 0.5)
    #plt.setp(stemlines, 'linestyle', 'dotted')
    
    if norm:
        markerline, stemlines, baseline = plt.stem(molecularmassTit,koncsorted/(konc[normVal]),markerfmt='o')
    else:
        markerline, stemlines, baseline = plt.stem(molecularmassTit,koncsorted,markerfmt='o')
    plt.setp(stemlines, 'color', plt.getp(markerline,'color'))
    plt.setp(stemlines, 'alpha', 0.7)
    plt.setp(markerline, 'alpha', 0.7)
    plt.setp(markerline, color='blue')
    plt.setp(stemlines, color='blue')
    plt.setp(stemlines, 'linestyle', 'dotted')
    
    #if norm:
     #   markerline, stemlines, baseline = plt.stem(molecularmass2,konc2sorted/(konc2[normVal]),markerfmt='o', label= "No Surfaces")
    #else:
     #   markerline, stemlines, baseline = plt.stem(molecularmass2,konc2sorted,markerfmt='o', label= "No Surfaces")
    #plt.setp(stemlines, 'color', plt.getp(markerline,'color'))
    #plt.setp(stemlines, 'alpha', 0.5)
    #plt.setp(markerline, 'alpha', 0.5)
    #plt.setp(markerline, color='red')
    #plt.setp(stemlines, color='red')
    #plt.setp(stemlines, 'linestyle', 'dotted')
    
    for label in labels:
        if val[0] < species[2][label[0]] < val[1]:
            if norm:
                plt.annotate(label[1], xy=(species[2][label[0]],  konc[label[0]]/(konc[normVal])), xycoords='data', xytext=(10, 20), textcoords='offset points', rotation=90, arrowprops=dict(arrowstyle="->", color='black'))
            else:
                plt.annotate(label[1], xy=(species[2][label[0]],  konc[label[0]]), xycoords='data', xytext=(10, 20), textcoords='offset points', rotation=90, arrowprops=dict(arrowstyle="->", color='black'))
    
    plt.xlabel("Molecular Mass (amu)")
    if norm:
        plt.ylabel(r"Intensity (au)")
    else:
        plt.ylabel(r"Concentration, n (cm$^{-3}$)")
    plt.xlim(val)
    if i==0:
        if norm:
            plt.ylim(1e-19,1e3)
        else:
            plt.ylim(1e-5,1e16)
    if i==1:
        plt.legend(loc='upper center')
    else:
        plt.legend()
    plt.yscale("log")
    plt.tight_layout()
    if norm:
        plt.savefig("TitAerFrTrial4"+str(i)+"Norm.pdf")
    else:
        plt.savefig("TitAerFrTrial4"+str(i)+".pdf")
    plt.show()