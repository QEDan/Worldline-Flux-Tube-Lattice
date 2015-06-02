import numpy as np
import matplotlib.pyplot as plt

f=open("/home/mazur/worldline/EA.F2l2p426.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I0=list()
xcm0=list()
eI0=list()
for row in varList:
	xcm0.append(eval(row[1]))
	I0.append(eval(row[2]))
	eI0.append(eval(row[3]))


f=open("/home/mazur/worldline/EA.F2l2p534.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I1=list()
xcm1=list()
eI1=list()
for row in varList:
	xcm1.append(eval(row[1]))
	I1.append(eval(row[2]))
	eI1.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F2l2p6.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I2=list()
xcm2=list()
eI2=list()
for row in varList:
	xcm2.append(eval(row[1]))
	I2.append(eval(row[2]))
	eI2.append(eval(row[3]))




plt.clf()
p0 =plt.errorbar(xcm0,I0,eI0,marker='.', color='k', label=r'$\lambda^2=0.426$',ms=8, linewidth=0.0, elinewidth=1)
p1 =plt.errorbar(xcm1,I1,eI1,marker='.', color='g', label=r'$\lambda^2=0.534$',ms=8, linewidth=0.0, elinewidth=1)
p2 =plt.errorbar(xcm2,I2,eI2,marker='.', color='r', label=r'$\lambda^2=0.6$', ms=8, linewidth=0.0, elinewidth=1)



#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Lagrange Density')
plt.legend()
ax=plt.gca()
#ax.set_ylim(-0.5,0.5)
ax.set_xlim(0.0,0.8)
plt.savefig('Lvr.F2.eps')
plt.savefig('Lvr.F2.png')
