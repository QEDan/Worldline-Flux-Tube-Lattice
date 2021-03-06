import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EAgvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l2=list()
EAclass=list()
EAferm=list()
eEAferm=list()
for row in varList:
	l2.append(eval(row[1])**2)
	EAclass.append(eval(row[2]))
	EAferm.append(eval(row[3]))
	eEAferm.append(eval(row[4]))

l=[sqrt(l2[i]) for i in range(len(l2))]

plt.clf()
plt.errorbar(l,[EAferm[i] + EAclass[i] for i in range(len(EAferm))],eEAferm,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.plot(l,EAclass,'k')
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Effective Action")
ax=plt.gca()
ax.set_ylim(-1000.0,2000)
#ax.set_xlim(0.0,0.8)
plt.savefig('EAgvsl.eps')
plt.savefig('EAgvsl.png')
