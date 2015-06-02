import numpy as np
import matplotlib.pyplot as plt

a = np.sqrt(8.0)

f=open("/home/mazur/worldline/EAscff.F1l2p16a2.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I=list()
xcm=list()
eI=list()
for row in varList:
	if row[0]=='Leffvsrho:':
		xcm.append(eval(row[1])/a)
		I.append(eval(row[2]))
		eI.append(eval(row[3]))

f=open("/home/mazur/worldline/EAscff.F1l2p25a2.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I2=list()
xcm2=list()
eI2=list()
for row in varList:
	if row[0]=='Leffvsrho:':
		xcm2.append(eval(row[1])/a)
		I2.append(eval(row[2]))
		eI2.append(eval(row[3]))

x=np.arange(0.0,40.0,0.1)

f=open("/home/mazur/worldline/EAscff.F1l2p36a2.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I3=list()
xcm3=list()
eI3=list()
for row in varList:
	if row[0]=='Leffvsrho:':
		xcm3.append(eval(row[1])/a)
		I3.append(eval(row[2]))
		eI3.append(eval(row[3]))


plt.clf()
plt.errorbar(xcm,I,eI,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$\lambda = 0.4a$')
plt.errorbar(xcm2,I2,eI2,marker='x', color='k', ms=4, linewidth=0.0, elinewidth=1, label=r'$\lambda = 0.5a$')
plt.errorbar(xcm3,I3,eI3,marker='^', color='k', ms=4, linewidth=0.0, elinewidth=1, label=r'$\lambda = 0.6a$')
#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}/a$")
plt.ylabel('Effective Action Density')
plt.legend(loc=0)
ax=plt.gca()
#ax.set_ylim(-0.1,0.1)
ax.set_xlim(0,0.30)
plt.savefig('EAscffpeek.png')
plt.savefig('EAscffpeek.eps')
