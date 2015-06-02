import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EApscvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l2=list()
EAclass=list()
EAferm=list()
eEAferm=list()
for row in varList:
	l2.append(eval(row[1])**2)
	EAclass.append(-eval(row[2]))
	EAferm.append(eval(row[3]))
	eEAferm.append(eval(row[4]))


data=zip(l2,EAclass,EAferm,eEAferm)
data.sort()
[l2,EAclass,EAferm,eEAferm]=zip(*data)
l=[sqrt(l2[i])for i in range(len(l2))]

plt.clf()
plt.errorbar(l,[EAferm[i]/EAclass[i] for i in range(len(EAferm))],[eEAferm[i]/EAclass[i] for i in range(len(EAferm))],marker='.', color='k',  ms=8, linewidth=1.0, elinewidth=1)
#plt.plot(l,EAclass,'k',linestyle='dashed')
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Effective Action")
ax=plt.gca()
#ax.set_yscale('log')

#ax.set_ylim(-1.0,0.4)
#ax.set_xlim(0.0,0.7)
plt.savefig('EApscvsl.eps')
plt.savefig('EApscvsl.png')
