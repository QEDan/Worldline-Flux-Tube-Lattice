import numpy as np
import matplotlib.pyplot as plt

m2=0.30282212;
e=0.30282212;
#m2=1.0;
#e=1.0;



f=open("/home/mazur/worldline/EAps.F1l2100a500.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I=list()
xcm=list()
eI=list()
for row in varList:
	xcm.append(eval(row[1]))
	I.append(eval(row[2]))
	eI.append(eval(row[3]))
x=np.arange(0.0,40.0,0.1)


plt.clf()
plt.errorbar(xcm,I,eI,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Action Density')
ax=plt.gca()
#ax.set_ylim(-1000,100)
#ax.set_xlim(0.0,0.8)
plt.savefig('EAps.F1l2100a500.eps')
plt.savefig('EAps.F1l2100a500.png')
