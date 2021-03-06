import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EA.Fp9l21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I1=list()
xcm1=list()
eI1=list()
for row in varList:
	xcm1.append(eval(row[1]))
	I1.append(eval(row[2]))
	eI1.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1l21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I2=list()
xcm2=list()
eI2=list()
for row in varList:
	xcm2.append(eval(row[1]))
	I2.append(eval(row[2]))
	eI2.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1p1l21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I3=list()
xcm3=list()
eI3=list()
for row in varList:
	xcm3.append(eval(row[1]))
	I3.append(eval(row[2]))
	eI3.append(eval(row[3]))




plt.clf()
p1 =plt.errorbar(xcm1,I1,eI1,marker='.', color='k', label=r'$\lambda^2=1, F=0.9$',ms=8, linewidth=0.0, elinewidth=1)
p2 =plt.errorbar(xcm2,I2,eI2,marker='x', color='k', label=r'$\lambda^2=1, F=1.0$',ms=8, linewidth=0.0, elinewidth=1)
p3 =plt.errorbar(xcm3,I3,eI3,marker='^', color='k', label=r'$\lambda^2=1, F=1.1$',ms=3, linewidth=0.0, elinewidth=1)



#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Lagrange Density')
plt.legend()
ax=plt.gca()
#ax.set_ylim(-0.5,0.5)
ax.set_xlim(0.0,0.8)
plt.savefig('Lvr.Fs.eps')
plt.savefig('Lvr.Fs.png')
