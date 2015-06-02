import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/ActRatsvsl.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l2=list()
ratio=list()
eratio=list()
for row in varList:
	l2.append(eval(row[1])**2)
	ratio.append(eval(row[2]))
	eratio.append(eval(row[3]))

plt.clf()
plt.errorbar(l2,ratio,eratio,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.xlabel(r"$\lambda^2$")
plt.ylabel(r"$\frac{\Gamma_{\rm ferm}}{\Gamma_{\rm class}}$")
ax=plt.gca()
ax.set_ylim(-0.00016,0.0)
ax.set_xlim(0.1,0.8)
plt.savefig('ActRatsvsl2.eps')
plt.savefig('ActRatsvsl2.png')
