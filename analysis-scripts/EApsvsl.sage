import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EApsvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l2=list()
EAferm=list()
eEAferm=list()
for row in varList:
	l2.append(eval(row[1])**2)
	EAferm.append(eval(row[3]))
	eEAferm.append(eval(row[4]))

plt.clf()
plt.errorbar(l2,EAferm,eEAferm,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.xlabel(r"$\lambda^2$")
plt.ylabel(r"$\Gamma^{(1)}_{\rm ferm}$")
ax=plt.gca()
ax.set_ylim(-3.0,0.4)
#ax.set_xlim(0.0,0.8)
plt.savefig('EApsvsl.eps')
plt.savefig('EApsvsl.png')
