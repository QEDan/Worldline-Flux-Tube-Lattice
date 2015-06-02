import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;

def B(rho2, l2, F):
	return F/e*l2/(l2+rho2)^2

def DELO(rho2, l2, F):
	var('T')
	f=numerical_integral(exp(-m2*T)/T**3*(l2*e*B(rho2,l2, F)*T/tanh(l2*e*B(rho2,l2, F)*T)
		-1.0-1.0/3.0*(l2*e*B(rho2,l2, F)*T)**2),(0,infinity))
	return f[0]

f=open("/home/mazur/worldline/EAscff.F1l2lmin2p25a210000.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I=list()
xcm=list()
eI=list()
for row in varList:
	if row[0]=='Leffvsrho:':
		xcm.append(eval(row[1]))
		I.append(eval(row[2]))
		eI.append(eval(row[3]))
x=np.arange(0.0,40.0,0.1)
DE=[DELO(xi**2,13.8, 1.0) for xi in x]

plt.clf()
plt.errorbar(xcm,I,eI,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
#plt.plot(x,DE,'k')
#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Action Density')
ax=plt.gca()
#ax.set_ylim(0,0.2)
ax.set_xlim(0.0,0.5)
plt.savefig('peek.png')
plt.savefig('peek.eps')
