import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;

def B(rho2,l2, F):
        return F/e*l2/(l2+rho2)^2


def EA(l2, F):
	return numerical_integral(lambda rho: 
	  rho*numerical_integral(lambda T: 
		exp(-m2*T)/T**3*(e*B(rho*rho,l2, F)*T/tanh(e*B(rho*rho,l2, F)*T)
		-1.0-1.0/3.0*(e*B(rho*rho,l2, F)*T)**2), 
          (0,infinity))[0],
	(0,infinity))[0]

	


f=open("/home/mazur/worldline/EAaf2vsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l=list()
EAferm=list()
eEAferm=list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	l.append(eval(row[1]))
	EAferm.append(-eval(row[3]))
	eEAferm.append(eval(row[4]))

f=open("/home/mazur/worldline/EAaf2vsLambda.F2.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l2=list()
EAferm2=list()
eEAferm2=list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	l2.append(eval(row[1]))
	EAferm2.append(-eval(row[3]))
	eEAferm2.append(eval(row[4]))

f=open("/home/mazur/worldline/EAaf2vsLambda.Fp5.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

lp5=list()
EAfermp5=list()
eEAfermp5=list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	lp5.append(eval(row[1]))
	EAfermp5.append(-eval(row[3]))
	eEAfermp5.append(eval(row[4]))

x=np.logspace(-1.4,1,14)
EffAct=[-EA(float(xi*xi),1.0) for xi in x]
EffAct2=[-EA(float(xi*xi),2.0) for xi in x]
EffActp5=[-EA(float(xi*xi),0.5) for xi in x]

plt.clf()
plt.errorbar(lp5,EAfermp5,eEAfermp5,marker='^', color='k', 
	ms=4, linewidth=0.0, elinewidth=1, label=r'$\mathcal{F} = 0.5$')
plt.errorbar(l,EAferm,eEAferm,marker='.', color='k', 
	ms=8, linewidth=0.0, elinewidth=1, label=r'$\mathcal{F} = 1$')
plt.errorbar(l2,EAferm2,eEAferm2,marker='x', color='k', 
	ms=4, linewidth=0.0, elinewidth=1, label=r'$\mathcal{F} = 2$')

plt.plot(x, EffActp5, 'k', linewidth=1.0, label=r'$\mathcal{F} = 0.5$')
plt.plot(x, EffAct, 'k--', linewidth=1.0, label=r'$\mathcal{F} = 1$')
plt.plot(x, EffAct2, 'k-.', linewidth=1.0, label=r'$\mathcal{F} = 2$')

plt.xlabel(r"$\lambda$")
plt.ylabel(r"$-\Gamma^{(1)}_{\rm ferm}$")
plt.legend()
ax=plt.gca()
plt.yscale('log')
plt.xscale('log')
#ax.set_ylim(1e-10,1e2)
#ax.set_xlim(0.0,0.8)
plt.savefig('EAaf2vsl.eps')
plt.savefig('EAaf2vsl.png')
