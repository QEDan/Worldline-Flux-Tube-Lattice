import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EAffvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l=list()
EAferm=list()
eEAferm=list()
for row in varList:
	l.append(eval(row[1])/sqrt(8.0))
	EAferm.append(eval(row[3]))
	eEAferm.append(eval(row[4]))

plt.clf()
plt.errorbar(l,EAferm,eEAferm,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.xlabel(r"$\lambda /a$")
plt.ylabel(r"$\Gamma^{(1)}_{\rm scal}$")
ax=plt.gca()
#ax.set_ylim(-0.03,0.04)
#ax.set_xlim(0.0,0.8)
plt.savefig('EAffvsl.eps')
plt.savefig('EAffvsl.png')
