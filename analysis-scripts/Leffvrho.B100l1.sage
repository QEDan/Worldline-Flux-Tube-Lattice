import numpy as np
import matplotlib.pyplot as plt

f=open("/home/mazur/worldline/EA.B100l1.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I=list()
xcm=list()
eI=list()
for row in varList:
	xcm.append(eval(row[1]))
	I.append(eval(row[2]))
	eI.append(eval(row[3]))

exact=0.0139688

plt.clf()
plt.errorbar(xcm,I,eI,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Lagrange Density')
plt.savefig('EffAct4.eps')
plt.savefig('EffAct4.png')
