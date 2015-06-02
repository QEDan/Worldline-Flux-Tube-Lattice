import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;
atol=1.0e-2;
rtol=1.0e-12;

call('grep EAvsLambda EAaf2.F1*.out > EAaf2vsLambda.dat', shell=True)

f=open("/home/mazur/worldline/EAaf2vsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

laf2=list()
EAaf2scal=list()
eEAaf2scal=list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	laf2.append(eval(row[1]))
	EAaf2scal.append(-eval(row[3]))
	eEAaf2scal.append(eval(row[4]))

call('grep EAvsLambda EAscsm*.out > EAscsmvsLambda.dat', shell=True)

f=open("/home/mazur/worldline/EAscsmvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

lscsm=list()
EAscsmscal=list()
eEAscsmscal=list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])) and eval(row[1]) >= 0.1:
	lscsm.append(eval(row[1]))
	EAscsmscal.append(-eval(row[3]))
	eEAscsmscal.append(eval(row[4]))

ratio = list()
eratio = list()
lratio = list()
for i in range(len(laf2)):
  for j in range(len(lscsm)):
    if abs(laf2[i] - lscsm[j]) < 1.0e-5:
	lratio.append(laf2[i])
	ratio.append(EAaf2scal[i]/EAscsmscal[j])
	eratio.append(EAaf2scal[i]/EAscsmscal[j]*np.sqrt(
	  (eEAaf2scal[i]/EAaf2scal[i])**2 + (eEAscsmscal[j]/EAscsmscal[j])**2)
	  -2.0*1.0*eEAaf2scal[i]*eEAscsmscal[j]/(EAaf2scal[i]*EAscsmscal[j]))
			 

def B(rho2,l2, F):
        return F/e*l2/(l2+rho2)**2


def EAferm(l2, F):
	return numerical_integral(lambda rho: 
	  rho*numerical_integral(lambda T: 
		exp(-m2*T)/T**3*(e*B(rho*rho,l2, F)*T/tanh(e*B(rho*rho,l2, F)*T)
		-1.0-1.0/3.0*(e*B(rho*rho,l2, F)*T)**2)/(4.0*pi), 
          (0,infinity), eps_abs=atol, eps_rel=rtol)[0],
	(0,infinity), eps_abs=atol, eps_rel=rtol)[0]

def EAscal(l2, F):
	return numerical_integral(lambda rho: 
	  rho*numerical_integral(lambda T: 
		-1*exp(-m2*T)/T**3*(e*B(rho*rho,l2, F)*T/sinh(e*B(rho*rho,l2, F)*T)
		-1.0+1.0/6.0*(e*B(rho*rho,l2, F)*T)**2)/(2.0*pi), 
          (0,100.0), eps_abs=atol, eps_rel=rtol)[0],
	(0,infinity), eps_abs=atol, eps_rel=rtol)[0]

#l = np.logspace(0,2,10)
#EffActferm=[abs(EAferm(float(li*li),1.0)) for li in l]
#EffActscal=[abs(EAscal(float(li*li),1.0)) for li in l]

f=open("/home/mazur/worldline/LODEferm.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l=list()
EffActferm = list()
for row in varList:
  l.append(eval(row[0]))
  EffActferm.append(-eval(row[1]))

f=open("/home/mazur/worldline/LODEscal.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l=list()
EffActscal = list()
for row in varList:
  l.append(eval(row[0]))
  EffActscal.append(-eval(row[1]))

plt.cla()
plt.plot(l,EffActferm,'k', label='Spinor QED')
plt.plot(l,EffActscal,'k--', label = 'Scalar QED')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$|\Gamma^{(1)}|$')
plt.xlabel(r'$\lambda$')
plt.legend()
plt.savefig('scalfermi.png')

EffActrat = list()
for i in range(len(EffActferm)):
  EffActrat.append(EffActferm[i]/EffActscal[i])

plt.cla()
plt.plot(l,EffActrat,'k',linewidth=2)
plt.errorbar(lratio,ratio, eratio, marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
#plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\Gamma^{(1)}_{\rm ferm}/\Gamma^{(1)}_{\rm scal}$', size='large')
plt.xlabel(r'$\lambda/\lambda_c$', size='large')
plt.savefig('scalfermiratio.png')
plt.savefig('scalfermiratio.eps')

call('scp scalfermiratio.eps mazur-desktop:~/Documents/Thesis/images/.', shell=True)
