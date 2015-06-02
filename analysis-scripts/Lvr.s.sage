import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;

def B(rho2,l2):
	return 1/(e*l2)*exp(-rho2/l2)

def DELO(rho2,l2):
	var('T')
	f=numerical_integral(exp(-m2*T)/T**3*(l2*e*B(rho2,l2)*T/sinh(l2*e*B(rho2,l2)*T)-1.0+1.0/6.0*(l2*e*B(rho2,l2)*T)**2),(0,infinity))
	return f[0]


f=open("/home/mazur/worldline/EAs.F1l210.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I1=list()
xcm1=list()
eI1=list()
for row in varList:
	xcm1.append(eval(row[1]))
	I1.append(eval(row[2]))
	eI1.append(eval(row[3]))

f=open("/home/mazur/worldline/EAs.F1l21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I2=list()
xcm2=list()
eI2=list()
for row in varList:
	xcm2.append(eval(row[1]))
	I2.append(eval(row[2]))
	eI2.append(eval(row[3]))

#DE=[DELO(xi**2,10.0)/DELO(0.0,10.0)*max(I1) for xi in np.sort(xcm1)]
DE=[max(I1)/B(0.0,10.0)**8*B(xi**2,10.0)**8 for xi in np.sort(xcm1)]


plt.clf()
p1 =plt.errorbar(xcm1,I1,eI1,marker='.', color='b', label=r'$\lambda^2=10$', ms=8, linewidth=0.0, elinewidth=1)
#p2 =plt.errorbar(xcm2,I2,eI2,marker='.', color='k', label=r'$\lambda^2=1$',ms=8, linewidth=0.0, elinewidth=1)
pe=plt.plot(np.sort(xcm1),DE)




#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Lagrange Density')
plt.legend()
ax=plt.gca()
#ax.set_ylim(-20.,200.)
ax.set_xlim(0.0,8.0)
plt.savefig('EAs.Lvr.eps')
plt.savefig('EAs.Lvr.png')
